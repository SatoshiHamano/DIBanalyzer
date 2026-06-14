#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Read-only database access helpers for DIB workbench tools.

All functions in this module are intended to run SELECT-only queries.  Keeping
them here lets the CLI inspector, profiler, and spectrum viewer share one DB
access path.
"""

from __future__ import print_function

import os
import sys


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)


def get_connection():
    from open_mysql_project import openproject

    return openproject()


def fetch_all_dict(cur):
    columns = [desc[0] for desc in cur.description]
    return [dict(zip(columns, row)) for row in cur.fetchall()]


def fetch_one_dict(cur):
    rows = fetch_all_dict(cur)
    if rows:
        return rows[0]
    return None


def objectdict_exact_lookup_query(query, limit):
    return (
        """
        select objectid
        from objectdict
        where registeredname = %s
        order by priority
        limit %s
        """,
        (query, limit),
    )


def object_rows_by_ids_query(object_ids):
    placeholders = ",".join(["%s"] * len(object_ids))
    return (
        """
        select
          o.objectid as object_id,
          o.objectname,
          o.type,
          o.objtype,
          o.sptype,
          o.ra,
          o.decli,
          o.E_BV as e_bv,
          group_concat(
            distinct od.registeredname order by od.priority separator ' | '
          ) as aliases
        from object o
        left join objectdict od on o.objectid = od.objectid
        where o.objectid in ({})
        group by
          o.objectid,
          o.objectname,
          o.type,
          o.objtype,
          o.sptype,
          o.ra,
          o.decli,
          o.E_BV
        order by o.objectid
        """.format(placeholders),
        tuple(object_ids),
    )


def object_search_like_query(query, limit):
    like = "%{}%".format(query)
    return (
        """
        select
          o.objectid as object_id,
          o.objectname,
          o.type,
          o.objtype,
          o.sptype,
          o.ra,
          o.decli,
          o.E_BV as e_bv,
          group_concat(
            distinct od.registeredname order by od.priority separator ' | '
          ) as aliases
        from object o
        left join objectdict od on o.objectid = od.objectid
        where o.objectname like %s or od.registeredname like %s
        group by
          o.objectid,
          o.objectname,
          o.type,
          o.objtype,
          o.sptype,
          o.ra,
          o.decli,
          o.E_BV
        order by o.objectid
        limit %s
        """,
        (like, like, limit),
    )


def object_detail_query(object_id):
    return (
        """
        select
          o.objectid as object_id,
          o.objectname,
          o.type,
          o.objtype,
          o.sptype,
          o.ra,
          o.decli,
          o.pmra,
          o.pmdec,
          o.epoch,
          o.Vmag,
          o.Jmag,
          o.Hmag,
          o.Kmag,
          o.E_BV as e_bv,
          o.radial_velocity,
          group_concat(
            distinct od.registeredname order by od.priority separator ' | '
          ) as aliases
        from object o
        left join objectdict od on o.objectid = od.objectid
        where o.objectid = %s
        group by
          o.objectid,
          o.objectname,
          o.type,
          o.objtype,
          o.sptype,
          o.ra,
          o.decli,
          o.pmra,
          o.pmdec,
          o.epoch,
          o.Vmag,
          o.Jmag,
          o.Hmag,
          o.Kmag,
          o.E_BV,
          o.radial_velocity
        """,
        (object_id,),
    )


def combine_list_query(object_id):
    return (
        """
        select
          combineID as combine_id,
          objectID as object_id,
          combinepath,
          mode,
          telluricflag,
          combineflag,
          combineNumber as combine_number
        from combinesummary
        where objectID = %s
        order by combineID
        """,
        (object_id,),
    )


def spectrum_orders_query(combine_id):
    return (
        """
        select
          combineID as combine_id,
          echelleorder as echelle_order,
          combinefilepath,
          lambdamin,
          lambdamax
        from combinedspectrum
        where combineID = %s
        order by echelleorder
        """,
        (combine_id,),
    )


def measurements_for_order_query(combine_id, order):
    return (
        """
        select
          m.measurementID as measurement_id,
          m.measurementNum as measurement_number,
          m.combineID as combine_id,
          m.DIBID as dib_id,
          d.wavelength_air,
          d.wavelength_vac,
          d.category,
          d.reference as dib_reference,
          m.echelleorder as echelle_order,
          m.DIBspecpath as dib_spec_path,
          m.DIBimgpath as dib_image_path,
          m.primaryflag,
          m.autonormalizeflag,
          m.automeasurementflag,
          m.multiorderflag,
          m.EW as ew,
          m.EWerr as ew_err,
          m.centerlam_air,
          m.centerlam_vac,
          m.helio_velocity,
          m.FWHM as fwhm,
          m.FWHMerr as fwhm_err,
          m.SNR as snr,
          m.integration_start,
          m.integration_end,
          m.peaklam_air,
          m.depth,
          m.comment
        from DIBmeasurement m
        join DIBlist d using(DIBID)
        where m.combineID = %s and m.echelleorder = %s
        order by d.wavelength_air, m.measurementID
        """,
        (combine_id, order),
    )


def dib_search_query(query, limit):
    like = "%{}%".format(query)
    return (
        """
        select
          DIBID as dib_id,
          wavelength_air,
          wavelength_vac,
          category,
          reference,
          FWHM as fwhm,
          wavenumber,
          comment
        from DIBlist
        where
          cast(DIBID as char) like %s
          or cast(wavelength_air as char) like %s
          or reference like %s
          or category like %s
        order by wavelength_air
        limit %s
        """,
        (like, like, like, like, limit),
    )


def infer_combine_id_query(object_id):
    return (
        """
        select combineID
        from combinesummary
        where objectID = %s
        order by combineID
        limit 1
        """,
        (object_id,),
    )


def object_search(cur, query, limit):
    sql, params = objectdict_exact_lookup_query(query, limit)
    cur.execute(sql, params)
    exact_ids = [row[0] for row in cur.fetchall()]
    if exact_ids:
        sql, params = object_rows_by_ids_query(exact_ids)
        cur.execute(sql, params)
        return fetch_all_dict(cur)

    sql, params = object_search_like_query(query, limit)
    cur.execute(sql, params)
    return fetch_all_dict(cur)


def object_detail(cur, object_id):
    sql, params = object_detail_query(object_id)
    cur.execute(sql, params)
    return fetch_one_dict(cur)


def combine_list(cur, object_id):
    sql, params = combine_list_query(object_id)
    cur.execute(sql, params)
    return fetch_all_dict(cur)


def spectrum_orders(cur, combine_id):
    sql, params = spectrum_orders_query(combine_id)
    cur.execute(sql, params)
    rows = fetch_all_dict(cur)
    for row in rows:
        path = row.get("combinefilepath")
        row["fits_exists"] = bool(path and os.path.exists(path))
    return rows


def measurements_for_order(cur, combine_id, order):
    sql, params = measurements_for_order_query(combine_id, order)
    cur.execute(sql, params)
    return fetch_all_dict(cur)


def dib_search(cur, query, limit):
    sql, params = dib_search_query(query, limit)
    cur.execute(sql, params)
    return fetch_all_dict(cur)


def infer_combine_id(cur, object_id):
    sql, params = infer_combine_id_query(object_id)
    cur.execute(sql, params)
    row = cur.fetchone()
    if row:
        return row[0]
    return None


def inspect_path(cur, object_id, combine_id, order):
    detail = object_detail(cur, object_id)
    combines = combine_list(cur, object_id)

    if combine_id is None and combines:
        combine_id = combines[0]["combine_id"]

    orders = []
    measurements = []
    selected_order = order

    if combine_id is not None:
        orders = spectrum_orders(cur, combine_id)
        if selected_order is None and orders:
            selected_order = orders[0]["echelle_order"]
        if selected_order is not None:
            measurements = measurements_for_order(cur, combine_id, selected_order)

    return {
        "object": detail,
        "combines": combines,
        "selected_combine_id": combine_id,
        "orders": orders,
        "selected_order": selected_order,
        "measurements": measurements,
        "summary": {
            "combine_count": len(combines),
            "order_count": len(orders),
            "measurement_count": len(measurements),
        },
    }
