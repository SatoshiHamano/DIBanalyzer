#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys


class Color:
    BLACK = "\033[30m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    BLUE = "\033[34m"
    END = "\033[0m"


def db_status(_args):
    from open_mysql_project import openproject

    conn, cur = openproject()
    cur.execute("select database(), version()")
    dbname, version = cur.fetchone()
    print("database: {}".format(dbname))
    print("mysql: {}".format(version))

    cur.execute(
        "select table_name, coalesce(table_rows, 0) "
        "from information_schema.tables "
        "where table_schema=database() order by table_name"
    )
    print("\ntables:")
    for table, rows in cur.fetchall():
        print("  {:24s} {:>8}".format(table, rows))

    cur.close()
    conn.close()
    return 0


def dib_status(args):
    from add_lineDIB_mysql import GetDIBList
    from open_mysql_project import openproject

    conn, cur = openproject()
    object_id = args.object_id

    cur.execute(
        "select objectname, type, sptype from object where objectid=%d;" % object_id
    )
    rows = cur.fetchall()
    if rows == []:
        print("Object ID {} is not found.".format(object_id))
        conn.close()
        return 1

    object_name, object_type, sp_type = rows[0]
    print("objectID: {}".format(object_id))
    print("object: {} ({}, {})".format(object_name, object_type, sp_type))

    dib_info = GetDIBList(args.wavelength_min, args.wavelength_max)
    dib_ids, wavelengths, _reference, categories, _fwhm, _wavenumber, _comment = dib_info

    print("\nDIB status:")
    missing = 0
    no_primary = 0
    multiple_primary = 0
    ok = 0

    for dib_id, wavelength, category in zip(dib_ids, wavelengths, categories):
        cur.execute(
            "SELECT x.measurementID, x.primaryflag "
            "from DIBmeasurement as x join combinesummary as y using(combineID) "
            "where y.objectID=%d and x.DIBID=%d;" % (object_id, dib_id)
        )
        rows = cur.fetchall()

        if rows == []:
            missing += 1
            if args.all:
                print(
                    Color.GREEN
                    + "DIB{:.1f} (ID:{}, {}): not measured".format(
                        wavelength, dib_id, category
                    )
                    + Color.END
                )
            continue

        measurement_ids = [i[0] for i in rows]
        primary_flags = [i[1] for i in rows]
        primary_count = sum(primary_flags)

        if primary_count == 1:
            ok += 1
            if args.all:
                primary_id = measurement_ids[primary_flags.index(1)]
                print(
                    Color.BLACK
                    + "DIB{:.1f} (ID:{}, {}): {} measured, primary={}".format(
                        wavelength, dib_id, category, len(measurement_ids), primary_id
                    )
                    + Color.END
                )
        elif primary_count > 1:
            multiple_primary += 1
            print(
                Color.RED
                + "DIB{:.1f} (ID:{}, {}): {} measured, multiple primary flags".format(
                    wavelength, dib_id, category, len(measurement_ids)
                )
                + Color.END
            )
        else:
            no_primary += 1
            print(
                Color.BLUE
                + "DIB{:.1f} (ID:{}, {}): {} measured, no primary flag".format(
                    wavelength, dib_id, category, len(measurement_ids)
                )
                + Color.END
            )

    print("\nsummary:")
    print("  ok: {}".format(ok))
    print("  not measured: {}".format(missing))
    print("  no primary: {}".format(no_primary))
    print("  multiple primary: {}".format(multiple_primary))

    cur.close()
    conn.close()
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        description="Small command-line entry point for DIBproject workflows."
    )
    subparsers = parser.add_subparsers(dest="command")

    db_parser = subparsers.add_parser("db", help="Database utilities")
    db_subparsers = db_parser.add_subparsers(dest="db_command")
    db_status_parser = db_subparsers.add_parser("status", help="Show DB status")
    db_status_parser.set_defaults(func=db_status)

    dib_parser = subparsers.add_parser("dib", help="DIB measurement utilities")
    dib_subparsers = dib_parser.add_subparsers(dest="dib_command")
    dib_status_parser = dib_subparsers.add_parser(
        "status", help="Show measurement status for one object"
    )
    dib_status_parser.add_argument("object_id", type=int)
    dib_status_parser.add_argument("--all", action="store_true", help="Show OK rows too")
    dib_status_parser.add_argument("--wavelength-min", type=float, default=9000.0)
    dib_status_parser.add_argument("--wavelength-max", type=float, default=13500.0)
    dib_status_parser.set_defaults(func=dib_status)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 2
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
