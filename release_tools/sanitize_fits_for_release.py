#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Create sanitized FITS copies for a data release.

The source FITS files are never modified.  This script reads
``release_inventory/candidate_release_inventory.csv``, removes private/internal header
keywords from each FITS header, and writes release copies under a separate
directory.
"""

import argparse
import csv
import hashlib
import os
import re

from astropy.io import fits


DEFAULT_REMOVE_KEYS = [
    "FITSFILE",
    "OBSERVER",
    "WODBPI",
    "WODBOBS",
    "WODBPROP",
]


def safe_name(value):
    value = str(value).strip().replace("*", "")
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"[^A-Za-z0-9_.+:-]+", "_", value)
    value = value.strip("_")
    return value or "UNKNOWN"


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sanitize_header(header, remove_keys, release_name, reference, doi):
    removed = []
    for key in remove_keys:
        if key in header:
            del header[key]
            removed.append(key)

    header["RELEASE"] = (release_name, "Data release package name")
    if reference:
        header["REFERENC"] = (reference[:68], "Reference publication")
    if doi:
        header["DOI"] = (doi[:68], "Data release DOI")
    header.add_history("Sanitized for public data release; source FITS unchanged.")
    if removed:
        header.add_history("Removed private/internal keys: {}".format(",".join(removed)))
    return removed


def output_path_for(row, output_dir):
    object_dir = safe_name(row.get("db_object_name") or row.get("index_star_name"))
    combine_id = safe_name(row.get("combine_id"))
    basename = os.path.basename(row["local_fits_path"])
    return os.path.join(output_dir, object_dir, combine_id, basename)


def sanitize_one(row, output_dir, remove_keys, release_name, reference, doi, overwrite):
    source = row["local_fits_path"]
    target = output_path_for(row, output_dir)
    if os.path.exists(target) and not overwrite:
        return target, [], "exists"

    os.makedirs(os.path.dirname(target), exist_ok=True)
    with fits.open(source) as hdul:
        removed_all = []
        for hdu in hdul:
            removed_all.extend(
                sanitize_header(hdu.header, remove_keys, release_name, reference, doi)
            )
        hdul.writeto(target, overwrite=overwrite, output_verify="silentfix")
    return target, sorted(set(removed_all)), "written"


def read_inventory(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def write_manifest(rows, manifest_path):
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
    fieldnames = [
        "release_candidate",
        "paper_2206_03131_role",
        "publication_candidate_tags",
        "publication_candidate_dois",
        "local_origin_tags",
        "local_origin_notes",
        "object_id",
        "db_object_name",
        "index_star_name",
        "combine_id",
        "echelle_order",
        "mode",
        "obs_date",
        "db_wavelength_min",
        "db_wavelength_max",
        "source_fits_name",
        "release_fits_path",
        "source_sha256",
        "release_sha256",
        "removed_header_keys",
        "status",
    ]
    with open(manifest_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def build_parser():
    parser = argparse.ArgumentParser(
        description="Create sanitized FITS copies for a DIB data release."
    )
    parser.add_argument(
        "--inventory",
        default="release_inventory/candidate_release_inventory.csv",
        help="Input inventory CSV.",
    )
    parser.add_argument(
        "--output-dir",
        default="release_package/spectra",
        help="Directory for sanitized FITS copies.",
    )
    parser.add_argument(
        "--manifest",
        default="release_package/MANIFEST_sanitized.csv",
        help="Output manifest for sanitized files.",
    )
    parser.add_argument(
        "--remove-key",
        action="append",
        default=[],
        help="Additional FITS header key to remove. Can be repeated.",
    )
    parser.add_argument(
        "--release-name",
        default="WINERED_DIB_PUBLISHED_DATA",
        help="Value written to RELEASE header keyword.",
    )
    parser.add_argument("--reference", default="", help="Publication reference.")
    parser.add_argument("--doi", default="", help="Dataset DOI, when available.")
    parser.add_argument("--limit", type=int, default=0, help="Only process first N rows.")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--no-checksum",
        action="store_true",
        help="Skip SHA256 checksums for faster test runs.",
    )
    return parser


def main():
    args = build_parser().parse_args()
    rows = read_inventory(args.inventory)
    if args.limit:
        rows = rows[: args.limit]

    remove_keys = DEFAULT_REMOVE_KEYS + args.remove_key
    manifest_rows = []

    for row in rows:
        source = row["local_fits_path"]
        target = output_path_for(row, args.output_dir)
        if args.dry_run:
            status = "dry-run"
            removed = []
        else:
            target, removed, status = sanitize_one(
                row,
                args.output_dir,
                remove_keys,
                args.release_name,
                args.reference,
                args.doi,
                args.overwrite,
            )

        source_sha = "" if args.no_checksum else sha256_file(source)
        release_sha = ""
        if not args.dry_run and os.path.exists(target) and not args.no_checksum:
            release_sha = sha256_file(target)

        manifest_rows.append(
            {
                "release_candidate": row.get("release_candidate", ""),
                "paper_2206_03131_role": row.get("paper_2206_03131_role", ""),
                "publication_candidate_tags": row.get("publication_candidate_tags", ""),
                "publication_candidate_dois": row.get("publication_candidate_dois", ""),
                "local_origin_tags": row.get("local_origin_tags", ""),
                "local_origin_notes": row.get("local_origin_notes", ""),
                "object_id": row.get("object_id", ""),
                "db_object_name": row.get("db_object_name", ""),
                "index_star_name": row.get("index_star_name", ""),
                "combine_id": row.get("combine_id", ""),
                "echelle_order": row.get("echelle_order", ""),
                "mode": row.get("mode", ""),
                "obs_date": row.get("obs_date", ""),
                "db_wavelength_min": row.get("db_wavelength_min", ""),
                "db_wavelength_max": row.get("db_wavelength_max", ""),
                "source_fits_name": os.path.basename(source),
                "release_fits_path": os.path.relpath(
                    os.path.abspath(target),
                    start=os.path.dirname(os.path.abspath(args.manifest)),
                ),
                "source_sha256": source_sha,
                "release_sha256": release_sha,
                "removed_header_keys": ",".join(removed),
                "status": status,
            }
        )

    if not args.dry_run:
        write_manifest(manifest_rows, args.manifest)

    print("input rows {}".format(len(rows)))
    print("dry_run {}".format(args.dry_run))
    print("output_dir {}".format(args.output_dir))
    print("manifest {}".format(args.manifest))
    print("remove_keys {}".format(",".join(remove_keys)))
    if not args.dry_run:
        print("written {}".format(sum(1 for row in manifest_rows if row["status"] == "written")))
        print("exists {}".format(sum(1 for row in manifest_rows if row["status"] == "exists")))


if __name__ == "__main__":
    main()
