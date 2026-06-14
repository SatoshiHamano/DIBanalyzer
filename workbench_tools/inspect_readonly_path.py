#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Inspect the read-only DB path needed by the first DIB workbench prototype.

This script intentionally runs SELECT queries only.  It is a small bridge
between the current MySQL schema and a future GUI/API layer.
"""

from __future__ import print_function

import argparse
import json
import os
import sys


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)


from workbench_tools.db_readonly import (  # noqa: E402
    dib_search,
    get_connection,
    inspect_path,
    object_search,
)


def print_json(payload):
    print(json.dumps(payload, indent=2, sort_keys=True, default=str))


def command_search_object(args):
    conn, cur = get_connection()
    try:
        print_json(object_search(cur, args.query, args.limit))
    finally:
        cur.close()
        conn.close()
    return 0


def command_search_dib(args):
    conn, cur = get_connection()
    try:
        print_json(dib_search(cur, args.query, args.limit))
    finally:
        cur.close()
        conn.close()
    return 0


def command_inspect(args):
    conn, cur = get_connection()
    try:
        print_json(inspect_path(cur, args.object_id, args.combine_id, args.order))
    finally:
        cur.close()
        conn.close()
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        description="Read-only DB path inspector for the DIB workbench prototype."
    )
    subparsers = parser.add_subparsers(dest="command")

    object_parser = subparsers.add_parser("object-search")
    object_parser.add_argument("query")
    object_parser.add_argument("--limit", type=int, default=20)
    object_parser.set_defaults(func=command_search_object)

    dib_parser = subparsers.add_parser("dib-search")
    dib_parser.add_argument("query")
    dib_parser.add_argument("--limit", type=int, default=20)
    dib_parser.set_defaults(func=command_search_dib)

    inspect_parser = subparsers.add_parser("inspect")
    inspect_parser.add_argument("object_id", type=int)
    inspect_parser.add_argument("--combine-id")
    inspect_parser.add_argument("--order", type=int)
    inspect_parser.set_defaults(func=command_inspect)

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
