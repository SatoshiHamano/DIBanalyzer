#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Profile read-only queries used by the DIB spectrum viewer.

This script runs SELECT and EXPLAIN only.  It is meant to identify obvious
query bottlenecks before adding any database indexes or migrations.
"""

from __future__ import print_function

import argparse
import json
import os
import sys
import time


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)


from workbench_tools.db_readonly import (  # noqa: E402
    combine_list_query,
    get_connection,
    infer_combine_id,
    measurements_for_order_query,
    object_detail_query,
    object_search_like_query,
    objectdict_exact_lookup_query,
    spectrum_orders_query,
)


def fetch_dicts(cur):
    columns = [desc[0] for desc in cur.description]
    return [dict(zip(columns, row)) for row in cur.fetchall()]


def run_query(cur, sql, params):
    cur.execute(sql, params)
    return cur.fetchall()


def time_query(cur, sql, params, repeat):
    times = []
    row_count = 0
    for _ in range(repeat):
        start = time.perf_counter()
        rows = run_query(cur, sql, params)
        times.append(time.perf_counter() - start)
        row_count = len(rows)
    return {
        "rows": row_count,
        "repeat": repeat,
        "min_ms": min(times) * 1000.0,
        "median_ms": sorted(times)[len(times) // 2] * 1000.0,
        "max_ms": max(times) * 1000.0,
    }


def explain_query(cur, sql, params):
    cur.execute("explain " + sql, params)
    return fetch_dicts(cur)


def table_indexes(cur, tables):
    placeholders = ",".join(["%s"] * len(tables))
    cur.execute(
        """
        select
          table_name,
          index_name,
          seq_in_index,
          column_name,
          non_unique,
          cardinality
        from information_schema.statistics
        where table_schema = database()
          and table_name in ({})
        order by table_name, index_name, seq_in_index
        """.format(placeholders),
        tuple(tables),
    )
    return fetch_dicts(cur)


def build_queries(object_id, object_query, combine_id, order):
    exact_sql, exact_params = objectdict_exact_lookup_query(object_query, 25)
    like_sql, like_params = object_search_like_query(object_query, 25)
    detail_sql, detail_params = object_detail_query(object_id)
    combine_sql, combine_params = combine_list_query(object_id)
    orders_sql, orders_params = spectrum_orders_query(combine_id)
    measurements_sql, measurements_params = measurements_for_order_query(combine_id, order)
    return [
        {
            "name": "objectdict_exact_lookup",
            "sql": exact_sql,
            "params": exact_params,
        },
        {
            "name": "object_search_like_fallback",
            "sql": like_sql,
            "params": like_params,
        },
        {
            "name": "object_detail",
            "sql": detail_sql,
            "params": detail_params,
        },
        {
            "name": "combine_list",
            "sql": combine_sql,
            "params": combine_params,
        },
        {
            "name": "spectrum_orders",
            "sql": orders_sql,
            "params": orders_params,
        },
        {
            "name": "measurements_for_order",
            "sql": measurements_sql,
            "params": measurements_params,
        },
    ]


def print_markdown(report):
    print("# Read-only query profile")
    print("")
    print("- object_id: {}".format(report["object_id"]))
    print("- object_query: {}".format(report["object_query"]))
    print("- combine_id: {}".format(report["combine_id"]))
    print("- order: {}".format(report["order"]))
    print("- repeat: {}".format(report["repeat"]))
    print("")

    print("## Timings")
    print("")
    print("| Query | Rows | Min ms | Median ms | Max ms |")
    print("| --- | ---: | ---: | ---: | ---: |")
    for item in report["queries"]:
        timing = item["timing"]
        print(
            "| `{}` | {} | {:.3f} | {:.3f} | {:.3f} |".format(
                item["name"],
                timing["rows"],
                timing["min_ms"],
                timing["median_ms"],
                timing["max_ms"],
            )
        )
    print("")

    print("## EXPLAIN summary")
    print("")
    print("| Query | Table | Type | Key | Rows | Extra |")
    print("| --- | --- | --- | --- | ---: | --- |")
    for item in report["queries"]:
        for row in item["explain"]:
            print(
                "| `{}` | `{}` | `{}` | `{}` | {} | {} |".format(
                    item["name"],
                    row.get("table"),
                    row.get("type"),
                    row.get("key"),
                    row.get("rows"),
                    row.get("Extra") or "",
                )
            )
    print("")

    print("## Existing indexes")
    print("")
    print("| Table | Index | Seq | Column | Non unique | Cardinality |")
    print("| --- | --- | ---: | --- | ---: | ---: |")
    for row in report["indexes"]:
        print(
            "| `{}` | `{}` | {} | `{}` | {} | {} |".format(
                row["table_name"],
                row["index_name"],
                row["seq_in_index"],
                row["column_name"],
                row["non_unique"],
                row["cardinality"],
            )
        )


def profile(args):
    conn, cur = get_connection()
    try:
        combine_id = args.combine_id or infer_combine_id(cur, args.object_id)
        if not combine_id:
            raise RuntimeError("No combineID found for object_id={}".format(args.object_id))

        report = {
            "object_id": args.object_id,
            "object_query": args.object_query,
            "combine_id": combine_id,
            "order": args.order,
            "repeat": args.repeat,
            "queries": [],
            "indexes": table_indexes(
                cur,
                [
                    "object",
                    "objectdict",
                    "combinesummary",
                    "combinedspectrum",
                    "DIBmeasurement",
                    "DIBlist",
                ],
            ),
        }

        for query in build_queries(args.object_id, args.object_query, combine_id, args.order):
            timing = time_query(cur, query["sql"], query["params"], args.repeat)
            explain = explain_query(cur, query["sql"], query["params"])
            report["queries"].append(
                {
                    "name": query["name"],
                    "timing": timing,
                    "explain": explain,
                }
            )

        if args.json:
            print(json.dumps(report, indent=2, sort_keys=True, default=str))
        else:
            print_markdown(report)
    finally:
        cur.close()
        conn.close()


def build_parser():
    parser = argparse.ArgumentParser(
        description="Profile SELECT-only queries used by the read-only workbench."
    )
    parser.add_argument("--object-id", type=int, default=59)
    parser.add_argument("--object-query", default="HD147889")
    parser.add_argument("--combine-id")
    parser.add_argument("--order", type=int, default=42)
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--json", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    profile(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
