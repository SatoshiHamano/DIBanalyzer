#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Compare existing telluric-corrected FITS with an Astropy recomputation.

This script is read-only for MySQL and input FITS files.  It writes validation
artifacts under ``telluric_validation/`` by default.
"""

from __future__ import print_function

import argparse
import csv
import glob
import os
import sys

import numpy as np
from astropy.io import fits
from astropy import units as u
from specutils import Spectrum1D
from specutils.manipulation import SplineInterpolatedResampler


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from open_mysql_project import openproject  # noqa: E402
from Spec1Dtools import openspecfits  # noqa: E402


def fetch_dict(cur):
    columns = [desc[0] for desc in cur.description]
    row = cur.fetchone()
    if row is None:
        return None
    return dict(zip(columns, row))


def fetch_telluric_row(cur, telluric_id=None, order=None):
    where = ["tr.telluricfilepath is not null", "tr.telluricfilepath != ''"]
    params = []
    if telluric_id:
        where.append("tr.telluricID = %s")
        params.append(telluric_id)
    if order is not None:
        where.append("tr.echelleorder = %s")
        params.append(order)

    sql = """
        select
          tr.telluricID,
          tr.echelleorder,
          tr.scale,
          tr.scaleerr,
          tr.shift,
          tr.shifterr,
          tr.lambdamin,
          tr.lambdamax,
          tr.frame,
          tr.telluricfilepath,
          tc.pipelineIDobj,
          tc.pipelineIDtel,
          tc.autoflag,
          tc.advanced,
          tc.pipelinever,
          tc.mode,
          tc.telluricNumber,
          tc.telluricPath
        from telluricresult tr
        join telluriccorrection tc using(telluricID)
        where {where}
        order by tr.telluricID desc, tr.echelleorder
        limit 1
    """.format(where=" and ".join(where))
    cur.execute(sql, tuple(params))
    return fetch_dict(cur)


def datareduction_path(cur, pipeline_id):
    cur.execute("select path from datareduction where pipelineID = %s", (pipeline_id,))
    row = cur.fetchone()
    if row:
        return row[0]
    return None


def infer_fsr_and_system(telluric_filepath):
    name = os.path.basename(telluric_filepath)
    parts = name.split("_")
    fsr = None
    system = None
    for part in parts:
        if part.startswith("fsr"):
            fsr = part
        if part in ("VAC", "AIR"):
            system = part
    if fsr is None or system is None:
        raise ValueError("Cannot infer fsr/system from {}".format(name))
    return fsr, system


def find_flux_fits(reduction_path, fsr, system, order, frame=None):
    if not reduction_path:
        return None
    patterns = []
    if frame:
        patterns.append(
            os.path.join(
                reduction_path,
                "*_{}".format(frame),
                "{}_flux".format(system),
                fsr,
                "*_m{}_*{}.fits".format(order, system),
            )
        )
    patterns.append(
        os.path.join(
            reduction_path,
            "*_sum",
            "{}_flux".format(system),
            fsr,
            "*_m{}_*{}.fits".format(order, system),
        )
    )
    pattern = os.path.join(
        reduction_path, "**", "*_m{}_*{}.fits".format(order, system)
    )
    patterns.append(pattern)
    for pattern in patterns:
        matches = sorted(glob.glob(pattern, recursive=True))
        if matches:
            return matches[0]
    return None


def recompute_telluric_astropy(inputsp, refsp, shift_pixel=0.0, scale=1.0, thres=0.001):
    inputspec = Spectrum1D.read(inputsp)
    refspec = Spectrum1D.read(refsp)
    inputfits = fits.open(inputsp)
    reffits = fits.open(refsp)
    try:
        aminput = float(inputfits[0].header["AIRMASS"])
        amref = float(reffits[0].header["AIRMASS"])
    finally:
        inputfits.close()
        reffits.close()

    dlamref = refspec.wavelength.value - np.roll(refspec.wavelength.value, 1)
    lamperpix = np.median(dlamref[1:])
    shiftedspec = Spectrum1D(
        spectral_axis=refspec.wavelength + lamperpix * shift_pixel * u.AA,
        flux=refspec.flux,
    )
    spline = SplineInterpolatedResampler(bin_edges="zero_fill")
    resampledspec = spline(shiftedspec, inputspec.wavelength)
    flux = np.array(resampledspec.flux.value, dtype=float)
    flux[flux < thres] = thres
    scaledspec = Spectrum1D(
        spectral_axis=resampledspec.wavelength,
        flux=flux * resampledspec.flux.unit,
    )
    scaledspec = Spectrum1D(
        spectral_axis=scaledspec.wavelength,
        flux=scaledspec.flux ** (aminput / amref * scale),
    )
    outputspec = inputspec / scaledspec
    return outputspec


def read_fits_xy(path):
    wavelength, flux, _crval, _cdelt, _crpix = openspecfits(path)
    return np.asarray(wavelength, dtype=float), np.asarray(flux, dtype=float)


def robust_stats(wavelength, existing_flux, recomputed_flux):
    finite = np.isfinite(existing_flux) & np.isfinite(recomputed_flux)
    if not finite.any():
        return {
            "n": 0,
            "median_diff": np.nan,
            "mad_diff": np.nan,
            "rms_diff": np.nan,
            "p95_abs_diff": np.nan,
            "max_abs_diff": np.nan,
            "median_frac_diff": np.nan,
        }
    diff = recomputed_flux[finite] - existing_flux[finite]
    denom = np.where(np.abs(existing_flux[finite]) > 0, existing_flux[finite], np.nan)
    frac = diff / denom
    return {
        "n": int(finite.sum()),
        "median_diff": float(np.nanmedian(diff)),
        "mad_diff": float(np.nanmedian(np.abs(diff - np.nanmedian(diff)))),
        "rms_diff": float(np.sqrt(np.nanmean(diff ** 2))),
        "p95_abs_diff": float(np.nanpercentile(np.abs(diff), 95)),
        "max_abs_diff": float(np.nanmax(np.abs(diff))),
        "median_frac_diff": float(np.nanmedian(frac)),
    }


def write_plot(output_png, wavelength, existing_flux, recomputed_flux, row, target_flux, ref_flux):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    diff = recomputed_flux - existing_flux
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    axes[0].plot(wavelength, existing_flux, lw=0.8, label="existing")
    axes[0].plot(wavelength, recomputed_flux, lw=0.8, alpha=0.75, label="recomputed")
    axes[0].set_ylabel("corrected flux")
    axes[0].legend(loc="best", fontsize=8)

    axes[1].plot(wavelength, diff, lw=0.8, color="tab:red")
    axes[1].axhline(0, color="0.3", lw=0.7)
    axes[1].set_ylabel("recomputed - existing")

    axes[2].plot(wavelength, target_flux, lw=0.7, label="target flux")
    axes[2].plot(wavelength, ref_flux, lw=0.7, label="ref flux", alpha=0.8)
    axes[2].set_ylabel("input flux")
    axes[2].set_xlabel("wavelength")
    axes[2].legend(loc="best", fontsize=8)

    fig.suptitle(
        "{} m{} scale={} shift={} A".format(
            row["telluricID"], row["echelleorder"], row["scale"], row["shift"]
        ),
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(output_png, dpi=150)
    plt.close(fig)


def validate(args):
    conn, cur = openproject()
    try:
        row = fetch_telluric_row(cur, args.telluric_id, args.order)
        if row is None:
            raise RuntimeError("No telluric result found")

        fsr, system = infer_fsr_and_system(row["telluricfilepath"])
        target_path = datareduction_path(cur, row["pipelineIDobj"])
        ref_path = datareduction_path(cur, row["pipelineIDtel"])
    finally:
        cur.close()
        conn.close()

    target_flux = find_flux_fits(target_path, fsr, system, row["echelleorder"], frame=row.get("frame"))
    ref_flux = find_flux_fits(ref_path, fsr, system, row["echelleorder"])
    existing = row["telluricfilepath"]

    missing = [
        name
        for name, path in [
            ("existing", existing),
            ("target_flux", target_flux),
            ("ref_flux", ref_flux),
        ]
        if not path or not os.path.exists(path)
    ]
    if missing:
        raise RuntimeError(
            "Missing paths: {}. existing={}, target_flux={}, ref_flux={}".format(
                ",".join(missing), existing, target_flux, ref_flux
            )
        )

    _target_wave, _target_flux_for_dw, _crval, target_dw, _crpix = openspecfits(target_flux)
    shift_pixel = float(row["shift"]) / float(target_dw)
    recomputed = recompute_telluric_astropy(
        target_flux,
        ref_flux,
        shift_pixel=shift_pixel,
        scale=float(row["scale"]),
        thres=args.threshold,
    )

    existing_wave, existing_flux = read_fits_xy(existing)
    recomputed_wave = recomputed.wavelength.value
    recomputed_flux = np.asarray(recomputed.flux.value, dtype=float)
    if len(existing_wave) != len(recomputed_wave) or np.nanmax(np.abs(existing_wave - recomputed_wave)) > 1e-4:
        recomputed_flux = np.interp(existing_wave, recomputed_wave, recomputed_flux)

    stats = robust_stats(existing_wave, existing_flux, recomputed_flux)
    os.makedirs(args.output_dir, exist_ok=True)
    safe_id = row["telluricID"].replace("/", "_")
    basename = "{}_m{}".format(safe_id, row["echelleorder"])
    output_png = os.path.join(args.output_dir, basename + "_astropy_compare.png")
    output_csv = os.path.join(args.output_dir, "summary.csv")

    target_wave, target_flux_arr = read_fits_xy(target_flux)
    ref_wave, ref_flux_arr = read_fits_xy(ref_flux)
    target_on_existing = np.interp(existing_wave, target_wave, target_flux_arr)
    ref_on_existing = np.interp(existing_wave, ref_wave, ref_flux_arr)
    if not args.no_plot:
        write_plot(output_png, existing_wave, existing_flux, recomputed_flux, row, target_on_existing, ref_on_existing)

    write_header = not os.path.exists(output_csv)
    with open(output_csv, "a", newline="") as handle:
        fieldnames = [
            "telluricID",
            "echelleorder",
            "system",
            "fsr",
            "scale",
            "shift_angstrom",
            "shift_pixel",
            "existing",
            "target_flux",
            "ref_flux",
            "plot",
        ] + sorted(stats.keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        payload = {
            "telluricID": row["telluricID"],
            "echelleorder": row["echelleorder"],
            "system": system,
            "fsr": fsr,
            "scale": row["scale"],
            "shift_angstrom": row["shift"],
            "shift_pixel": shift_pixel,
            "existing": existing,
            "target_flux": target_flux,
            "ref_flux": ref_flux,
            "plot": output_png if not args.no_plot else "",
        }
        payload.update(stats)
        writer.writerow(payload)

    print("telluricID: {}".format(row["telluricID"]))
    print("order: {}".format(row["echelleorder"]))
    print("system/fsr: {} {}".format(system, fsr))
    print("existing: {}".format(existing))
    print("target_flux: {}".format(target_flux))
    print("ref_flux: {}".format(ref_flux))
    print("shift: {:.6f} A = {:.6f} pix".format(float(row["shift"]), shift_pixel))
    print("scale: {}".format(row["scale"]))
    for key in sorted(stats.keys()):
        print("{}: {}".format(key, stats[key]))
    if not args.no_plot:
        print("plot: {}".format(output_png))
    print("summary: {}".format(output_csv))


def build_parser():
    parser = argparse.ArgumentParser(description="Validate Astropy telluric recomputation against existing output.")
    parser.add_argument("--telluric-id")
    parser.add_argument("--order", type=int, default=42)
    parser.add_argument("--threshold", type=float, default=0.001)
    parser.add_argument("--output-dir", default="telluric_validation")
    parser.add_argument("--no-plot", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    validate(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
