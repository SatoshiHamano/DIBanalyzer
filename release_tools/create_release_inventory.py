#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Create a local inventory for the candidate DIB spectra release set.

This reads ``order_spec_index.xlsx`` without external Excel dependencies,
matches the combine IDs against the MySQL database, and writes a CSV inventory.
It does not copy spectra or modify the database.
"""

import argparse
import csv
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from open_mysql_project import openproject


XLSX_NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

PAPER_2206_03131_REDDENED_OBJECT_IDS = {
    50,
    75,
    23,
    22,
    24,
    33,
    43,
    30,
    123,
    31,
    42,
    124,
    90,
    246,
    32,
    27,
    63,
    45,
    59,
    41,
    66,
    40,
    39,
    14,
    10,
    12,
    13,
    15,
    9,
    132,
    144,
}

PAPER_2206_03131_REFERENCE_OBJECT_IDS = {138}

PAPER_2016_CYGOB2_OBJECT_IDS = {
    9,   # Cyg OB2 No.12
    10,  # Cyg OB2 No.10
    11,  # Cyg OB2 No.11
    12,  # Cyg OB2 No.3
    13,  # Cyg OB2 No.5 / BD+40 4220
    14,  # Cyg OB2 No.8A
    15,  # Cyg OB2 No.9
}

PAPER_2016_CYGOB2_DOI = "10.3847/0004-637X/821/1/42"

PAPER_2015_NIR_DIB_OBJECT_IDS = {
    5,    # HD 202850
    9,    # Cyg OB2 No.12
    55,   # HD 12953
    57,   # HD 14489
    68,   # HD 190603
    75,   # HD 20041
    85,   # HD 21291
    86,   # HD 21389
    90,   # HD 223385
    94,   # HD 23180
    97,   # HD 24398
    100,  # HD 24912
    101,  # HD 25204
    106,  # HD 2905
    108,  # HD 30614
    110,  # HD 36371
    111,  # HD 36486
    112,  # HD 36822
    117,  # HD 37043
    119,  # HD 37128
    120,  # HD 37742
    122,  # HD 38771
    123,  # HD 41117
    124,  # HD 43384
    128,  # HD 50064
}

PAPER_2015_NIR_DIB_REFERENCE_OBJECT_IDS = {138}
PAPER_2015_NIR_DIB_DOI = "10.1088/0004-637X/800/2/137"

NTT17C_ANALYSIS_CANDIDATE_OBJECT_IDS = {
    8,
    25,
    26,
    28,
    29,
    34,
    35,
    37,
    44,
    52,
    60,
    65,
    133,
    135,
}


def _column_number(cell_ref):
    match = re.match(r"([A-Z]+)", cell_ref)
    if match is None:
        return 0
    number = 0
    for char in match.group(1):
        number = number * 26 + ord(char) - 64
    return number - 1


def _cell_value(cell, shared_strings):
    cell_type = cell.attrib.get("t")
    value = cell.find("a:v", XLSX_NS)

    if cell_type == "s" and value is not None:
        return shared_strings[int(value.text)]
    if cell_type == "inlineStr":
        return "".join(text.text or "" for text in cell.findall(".//a:t", XLSX_NS))
    if value is not None:
        return value.text or ""
    return ""


def read_order_index(path):
    with zipfile.ZipFile(path) as archive:
        shared_strings = []
        try:
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for string_item in root.findall("a:si", XLSX_NS):
                shared_strings.append(
                    "".join(text.text or "" for text in string_item.findall(".//a:t", XLSX_NS))
                )
        except KeyError:
            pass

        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))

    rows = []
    for row in sheet.findall(".//a:row", XLSX_NS):
        values = {}
        for cell in row.findall("a:c", XLSX_NS):
            values[_column_number(cell.attrib.get("r", "A"))] = _cell_value(
                cell, shared_strings
            )
        if values:
            rows.append([values.get(index, "") for index in range(8)])

    header = rows[0][1:]
    header_index = {name: index for index, name in enumerate(header)}
    records = []
    for raw_row in rows[1:]:
        row = raw_row[1:]
        spec_path = row[header_index["spec_path"]]
        parts = spec_path.split("/")
        combine_id = parts[1] if len(parts) >= 2 else ""
        records.append(
            {
                "index_setting": row[header_index["setting"]],
                "echelle_order": int(float(row[header_index["order"]])),
                "index_star_name": row[header_index["star_name"]],
                "obs_date": row[header_index["obs_date"]],
                "index_spec_path": spec_path,
                "index_object_dir": parts[0] if len(parts) >= 1 else "",
                "combine_id": combine_id,
                "index_x_min": row[header_index["x_min"]],
                "index_x_max": row[header_index["x_max"]],
            }
        )
    return records


def fetch_db_rows(combine_ids):
    conn, cur = openproject()
    placeholders = ",".join(["%s"] * len(combine_ids))
    query = (
        "select "
        "cs.combineID, cs.echelleorder, cs.combinefilepath, cs.lambdamin, cs.lambdamax, "
        "c.objectID, o.objectname, o.sptype, o.ra, o.decli, o.E_BV, "
        "c.combinepath, c.mode, c.telluricflag, c.combineflag, "
        "group_concat(od.registeredname order by od.priority separator ' | ') "
        "from combinedspectrum cs "
        "join combinesummary c on cs.combineID=c.combineID "
        "left join object o on c.objectID=o.objectid "
        "left join objectdict od on c.objectID=od.objectid "
        "where cs.combineID in ({})"
        "group by cs.combineID, cs.echelleorder, cs.combinefilepath, cs.lambdamin, cs.lambdamax, "
        "c.objectID, o.objectname, o.sptype, o.ra, o.decli, o.E_BV, "
        "c.combinepath, c.mode, c.telluricflag, c.combineflag "
    ).format(placeholders)
    cur.execute(query, combine_ids)
    rows = cur.fetchall()
    cur.close()
    conn.close()

    db_rows = {}
    for row in rows:
        key = (row[0], int(row[1]))
        db_rows[key] = {
            "db_combine_id": row[0],
            "db_echelle_order": int(row[1]),
            "local_fits_path": row[2],
            "db_wavelength_min": row[3],
            "db_wavelength_max": row[4],
            "object_id": row[5],
            "db_object_name": row[6],
            "db_sptype": row[7],
            "db_ra": row[8],
            "db_dec": row[9],
            "db_e_bv": row[10],
            "combine_path": row[11],
            "mode": row[12],
            "telluricflag": row[13],
            "combineflag": row[14],
            "db_object_aliases": row[15],
        }
    return db_rows


def build_inventory(index_path):
    order_rows = read_order_index(index_path)
    combine_ids = sorted({row["combine_id"] for row in order_rows})
    db_rows = fetch_db_rows(combine_ids)

    inventory = []
    for row in order_rows:
        key = (row["combine_id"], row["echelle_order"])
        db_row = db_rows.get(key, {})
        local_path = db_row.get("local_fits_path", "")
        file_exists = bool(local_path and os.path.exists(local_path))
        file_size_bytes = os.path.getsize(local_path) if file_exists else ""
        object_id = db_row.get("object_id")
        if object_id in PAPER_2206_03131_REDDENED_OBJECT_IDS:
            paper_role = "paper_2206_03131_reddened_target"
        elif object_id in PAPER_2206_03131_REFERENCE_OBJECT_IDS:
            paper_role = "paper_2206_03131_reference_star"
        else:
            paper_role = "extra_candidate_release_object"

        publication_tags = []
        publication_dois = []
        if object_id in PAPER_2016_CYGOB2_OBJECT_IDS:
            publication_tags.append("Hamano2016_CygOB2_DIB")
            publication_dois.append(PAPER_2016_CYGOB2_DOI)
        if object_id in PAPER_2015_NIR_DIB_OBJECT_IDS:
            publication_tags.append("Hamano2015_NIR_DIB_0p91_1p32um")
            publication_dois.append(PAPER_2015_NIR_DIB_DOI)
        if object_id in PAPER_2015_NIR_DIB_REFERENCE_OBJECT_IDS:
            publication_tags.append("Hamano2015_NIR_DIB_reference_star")
            publication_dois.append(PAPER_2015_NIR_DIB_DOI)

        local_origin_tags = []
        local_origin_notes = []
        if object_id in NTT17C_ANALYSIS_CANDIDATE_OBJECT_IDS:
            local_origin_tags.append("NTT17c_DIBanalysis_candidate")
            local_origin_notes.append(
                "Found in 2017-07 NTT17c quality/telluric/DIBanalysis traces; not yet tied to a publication."
            )

        inventory.append(
            {
                **row,
                **db_row,
                "paper_2206_03131_role": paper_role,
                "publication_candidate_tags": ";".join(publication_tags),
                "publication_candidate_dois": ";".join(publication_dois),
                "local_origin_tags": ";".join(local_origin_tags),
                "local_origin_notes": " ".join(local_origin_notes),
                "file_exists": file_exists,
                "file_size_bytes": file_size_bytes,
                "release_candidate": file_exists and bool(db_row),
            }
        )
    return inventory


def write_inventory(rows, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fieldnames = [
        "release_candidate",
        "paper_2206_03131_role",
        "publication_candidate_tags",
        "publication_candidate_dois",
        "local_origin_tags",
        "local_origin_notes",
        "object_id",
        "db_object_name",
        "db_object_aliases",
        "db_sptype",
        "db_ra",
        "db_dec",
        "db_e_bv",
        "index_star_name",
        "combine_id",
        "echelle_order",
        "mode",
        "telluricflag",
        "combineflag",
        "obs_date",
        "db_wavelength_min",
        "db_wavelength_max",
        "index_x_min",
        "index_x_max",
        "local_fits_path",
        "file_exists",
        "file_size_bytes",
        "index_spec_path",
        "index_object_dir",
        "combine_path",
        "index_setting",
    ]
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", default="order_spec_index.xlsx")
    parser.add_argument(
        "--output",
        default="release_inventory/candidate_release_inventory.csv",
        help="Output CSV path.",
    )
    args = parser.parse_args()

    inventory = build_inventory(args.index)
    write_inventory(inventory, args.output)

    print("wrote {}".format(args.output))
    print("rows {}".format(len(inventory)))
    print("release candidates {}".format(sum(1 for row in inventory if row["release_candidate"])))
    print("missing files {}".format(sum(1 for row in inventory if not row["file_exists"])))
    print("objects {}".format(len({row["object_id"] for row in inventory})))
    print("combines {}".format(len({row["combine_id"] for row in inventory})))


if __name__ == "__main__":
    main()
