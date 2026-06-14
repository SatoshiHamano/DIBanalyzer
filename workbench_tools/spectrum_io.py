#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Spectrum FITS reading helpers for DIB workbench tools."""

from __future__ import print_function

import os

import numpy as np
from astropy.io import fits


def downsample_xy(x, y, max_points):
    """Downsample a spectrum while preserving min/max flux in each bin."""
    if max_points <= 0 or len(x) <= max_points:
        return x, y

    bins = int(max_points / 2)
    if bins < 1:
        bins = 1
    edges = np.linspace(0, len(x), bins + 1, dtype=int)
    xs = []
    ys = []

    for start, end in zip(edges[:-1], edges[1:]):
        if end <= start:
            continue
        chunk_y = y[start:end]
        finite = np.isfinite(chunk_y)
        if not finite.any():
            continue
        chunk_x = x[start:end]
        local_y = chunk_y[finite]
        local_x = chunk_x[finite]
        min_pos = int(np.argmin(local_y))
        max_pos = int(np.argmax(local_y))
        pair = sorted(
            [
                (local_x[min_pos], local_y[min_pos]),
                (local_x[max_pos], local_y[max_pos]),
            ],
            key=lambda item: item[0],
        )
        for px, py in pair:
            xs.append(float(px))
            ys.append(float(py))

    return np.array(xs), np.array(ys)


def wavelength_from_header(header, size):
    """Build a linear wavelength array from common 1D FITS WCS keywords."""
    crval = header.get("CRVAL1")
    crpix = header.get("CRPIX1", 1.0)
    cdelt = header.get("CDELT1", header.get("CD1_1"))
    if crval is None or cdelt is None:
        return np.arange(size, dtype=float)
    pixel = np.arange(size, dtype=float) + 1.0
    return float(crval) + (pixel - float(crpix)) * float(cdelt)


def read_primary_1d(path):
    """Read the primary HDU as a flattened 1D spectrum and wavelength array."""
    if not path or not os.path.exists(path):
        raise ValueError("FITS file is not found: {}".format(path))

    with fits.open(path, memmap=False) as hdul:
        if not hdul or hdul[0].data is None:
            raise ValueError("Primary HDU has no data: {}".format(path))
        data = np.asarray(hdul[0].data, dtype=float).reshape(-1)
        wavelength = wavelength_from_header(hdul[0].header, data.size)

    return wavelength, data


def filter_spectrum(wavelength, flux, wavelength_min=None, wavelength_max=None):
    mask = np.isfinite(wavelength) & np.isfinite(flux)
    if wavelength_min is not None:
        mask &= wavelength >= wavelength_min
    if wavelength_max is not None:
        mask &= wavelength <= wavelength_max
    return wavelength[mask], flux[mask]


def spectrum_payload(path, max_points, wavelength_min=None, wavelength_max=None):
    """Read a FITS spectrum and return JSON-serializable arrays and metadata."""
    wavelength, flux = read_primary_1d(path)
    wavelength, flux = filter_spectrum(
        wavelength,
        flux,
        wavelength_min=wavelength_min,
        wavelength_max=wavelength_max,
    )

    if wavelength.size == 0:
        return {
            "wavelength": [],
            "flux": [],
            "points": 0,
            "full_points": 0,
            "downsampled": False,
        }

    full_points = int(wavelength.size)
    wavelength, flux = downsample_xy(wavelength, flux, max_points)

    return {
        "wavelength": [float(v) for v in wavelength],
        "flux": [float(v) for v in flux],
        "points": int(len(wavelength)),
        "full_points": full_points,
        "downsampled": full_points != len(wavelength),
    }


def robust_flux_limits(flux, low=0.01, high=0.99):
    """Return percentile-based y-limits for display code or diagnostics."""
    values = np.asarray(flux, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0, 1.0
    ymin, ymax = np.quantile(values, [low, high])
    if ymin == ymax:
        ymin -= 0.05
        ymax += 0.05
    pad = (ymax - ymin) * 0.08
    return float(ymin - pad), float(ymax + pad)
