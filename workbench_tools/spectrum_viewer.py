#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Small read-only spectrum viewer for DIBproject.

The server uses only the Python standard library for HTTP, plus the existing
project DB connection and astropy/numpy for FITS reading.  All DB operations are
SELECT-only.
"""

from __future__ import print_function

import argparse
import json
import os
import sys
from http.server import BaseHTTPRequestHandler, HTTPServer
from socketserver import ThreadingMixIn
from urllib.parse import parse_qs, urlparse

import numpy as np


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from workbench_tools.db_readonly import (  # noqa: E402
    combine_list,
    get_connection,
    measurements_for_order,
    object_detail,
    object_search,
    spectrum_orders,
)
from workbench_tools.review_notes import (  # noqa: E402
    STATUSES,
    selected_notes,
    upsert_note,
)
from workbench_tools.spectrum_io import spectrum_payload  # noqa: E402


def json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    return str(value)


def parse_int(value, default=None):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def parse_float(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def get_first(values, key, default=None):
    value = values.get(key, [default])[0]
    if value == "":
        return default
    return value


def response_json(handler, payload, status=200):
    body = json.dumps(payload, default=json_default, sort_keys=True).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def response_html(handler, body):
    data = body.encode("utf-8")
    handler.send_response(200)
    handler.send_header("Content-Type", "text/html; charset=utf-8")
    handler.send_header("Content-Length", str(len(data)))
    handler.end_headers()
    handler.wfile.write(data)


class ThreadedHTTPServer(ThreadingMixIn, HTTPServer):
    daemon_threads = True


class ViewerHandler(BaseHTTPRequestHandler):
    server_version = "DIBSpectrumViewer/0.1"

    def log_message(self, fmt, *args):
        sys.stderr.write("%s - - [%s] %s\n" % (self.address_string(), self.log_date_time_string(), fmt % args))

    def do_GET(self):
        parsed = urlparse(self.path)
        query = parse_qs(parsed.query)

        try:
            if parsed.path == "/":
                response_html(self, HTML)
            elif parsed.path == "/api/objects":
                self.handle_objects(query)
            elif parsed.path == "/api/object":
                self.handle_object(query)
            elif parsed.path == "/api/orders":
                self.handle_orders(query)
            elif parsed.path == "/api/spectrum":
                self.handle_spectrum(query)
            elif parsed.path == "/api/compare":
                self.handle_compare(query)
            elif parsed.path == "/api/compare-combines":
                self.handle_compare_combines(query)
            elif parsed.path == "/api/review-notes":
                self.handle_review_notes(query)
            elif parsed.path == "/api/health":
                response_json(self, {"ok": True, "mode": "read-only"})
            else:
                response_json(self, {"error": "not found"}, status=404)
        except Exception as exc:
            response_json(self, {"error": str(exc)}, status=500)

    def do_POST(self):
        parsed = urlparse(self.path)
        try:
            if parsed.path == "/api/review-note":
                self.handle_save_review_note()
            else:
                response_json(self, {"error": "not found"}, status=404)
        except Exception as exc:
            response_json(self, {"error": str(exc)}, status=500)

    def handle_objects(self, query):
        q = get_first(query, "q", "")
        limit = parse_int(get_first(query, "limit", "20"), 20)
        conn, cur = get_connection()
        try:
            response_json(self, object_search(cur, q, limit))
        finally:
            cur.close()
            conn.close()

    def handle_object(self, query):
        object_id = parse_int(get_first(query, "object_id"))
        if object_id is None:
            response_json(self, {"error": "object_id is required"}, status=400)
            return

        conn, cur = get_connection()
        try:
            detail = object_detail(cur, object_id)
            combines = combine_list(cur, object_id)
            response_json(self, {"object": detail, "combines": combines})
        finally:
            cur.close()
            conn.close()

    def handle_orders(self, query):
        combine_id = get_first(query, "combine_id")
        if not combine_id:
            response_json(self, {"error": "combine_id is required"}, status=400)
            return

        conn, cur = get_connection()
        try:
            response_json(self, spectrum_orders(cur, combine_id))
        finally:
            cur.close()
            conn.close()

    def handle_spectrum(self, query):
        combine_id = get_first(query, "combine_id")
        order = parse_int(get_first(query, "order"))
        max_points = parse_int(get_first(query, "max_points", "2500"), 2500)
        wave_min = parse_float(get_first(query, "wave_min"))
        wave_max = parse_float(get_first(query, "wave_max"))

        if not combine_id or order is None:
            response_json(self, {"error": "combine_id and order are required"}, status=400)
            return

        conn, cur = get_connection()
        try:
            orders = spectrum_orders(cur, combine_id)
            selected = None
            for item in orders:
                if int(item["echelle_order"]) == int(order):
                    selected = item
                    break
            if selected is None:
                response_json(self, {"error": "order is not found for combine_id"}, status=404)
                return

            measurements = measurements_for_order(cur, combine_id, order)
            spectrum = spectrum_payload(
                selected["combinefilepath"],
                max_points=max_points,
                wavelength_min=wave_min,
                wavelength_max=wave_max,
            )
            response_json(
                self,
                {
                    "combine_id": combine_id,
                    "order": order,
                    "spectrum_file": selected,
                    "spectrum": spectrum,
                    "measurements": measurements,
                },
            )
        finally:
            cur.close()
            conn.close()

    def handle_compare(self, query):
        object_ids_raw = get_first(query, "object_ids", "")
        wavelength = parse_float(get_first(query, "wavelength"))
        window = parse_float(get_first(query, "window", "30.0"), 30.0)
        max_points = parse_int(get_first(query, "max_points", "900"), 900)

        object_ids = [
            parse_int(item)
            for item in object_ids_raw.replace(" ", ",").split(",")
            if item.strip()
        ]
        object_ids = [item for item in object_ids if item is not None]

        if not object_ids or wavelength is None:
            response_json(self, {"error": "object_ids and wavelength are required"}, status=400)
            return

        wave_min = wavelength - window
        wave_max = wavelength + window
        rows = []

        conn, cur = get_connection()
        try:
            for object_id in object_ids:
                detail = object_detail(cur, object_id)
                if not detail:
                    rows.append({"object_id": object_id, "error": "object is not found"})
                    continue

                combines = combine_list(cur, object_id)
                if not combines:
                    rows.append({"object": detail, "error": "combine is not found"})
                    continue

                combine = combines[0]
                orders = spectrum_orders(cur, combine["combine_id"])
                selected = None
                for order in orders:
                    if float(order["lambdamin"]) <= wavelength <= float(order["lambdamax"]):
                        selected = order
                        break
                if selected is None:
                    rows.append({
                        "object": detail,
                        "combine": combine,
                        "error": "no order covers wavelength",
                    })
                    continue

                measurements = measurements_for_order(
                    cur,
                    combine["combine_id"],
                    selected["echelle_order"],
                )
                spectrum = spectrum_payload(
                    selected["combinefilepath"],
                    max_points=max_points,
                    wavelength_min=wave_min,
                    wavelength_max=wave_max,
                )
                rows.append(
                    {
                        "object": detail,
                        "combine": combine,
                        "order": selected,
                        "spectrum": spectrum,
                        "measurements": measurements,
                    }
                )

            response_json(
                self,
                {
                    "object_ids": object_ids,
                    "wavelength": wavelength,
                    "window": window,
                    "wave_min": wave_min,
                    "wave_max": wave_max,
                    "rows": rows,
                },
            )
        finally:
            cur.close()
            conn.close()

    def handle_compare_combines(self, query):
        object_id = parse_int(get_first(query, "object_id"))
        wavelength = parse_float(get_first(query, "wavelength"))
        window = parse_float(get_first(query, "window", "30.0"), 30.0)
        max_points = parse_int(get_first(query, "max_points", "900"), 900)

        if object_id is None or wavelength is None:
            response_json(self, {"error": "object_id and wavelength are required"}, status=400)
            return

        wave_min = wavelength - window
        wave_max = wavelength + window
        rows = []

        conn, cur = get_connection()
        try:
            detail = object_detail(cur, object_id)
            if not detail:
                response_json(self, {"error": "object is not found"}, status=404)
                return

            combines = combine_list(cur, object_id)
            for combine in combines:
                orders = spectrum_orders(cur, combine["combine_id"])
                selected = None
                for order in orders:
                    if float(order["lambdamin"]) <= wavelength <= float(order["lambdamax"]):
                        selected = order
                        break
                if selected is None:
                    rows.append({
                        "object": detail,
                        "combine": combine,
                        "error": "no order covers wavelength",
                    })
                    continue

                measurements = measurements_for_order(
                    cur,
                    combine["combine_id"],
                    selected["echelle_order"],
                )
                spectrum = spectrum_payload(
                    selected["combinefilepath"],
                    max_points=max_points,
                    wavelength_min=wave_min,
                    wavelength_max=wave_max,
                )
                rows.append(
                    {
                        "object": detail,
                        "combine": combine,
                        "order": selected,
                        "spectrum": spectrum,
                        "measurements": measurements,
                    }
                )

            response_json(
                self,
                {
                    "object_id": object_id,
                    "object": detail,
                    "wavelength": wavelength,
                    "window": window,
                    "wave_min": wave_min,
                    "wave_max": wave_max,
                    "rows": rows,
                    "comparison_type": "combines",
                },
            )
        finally:
            cur.close()
            conn.close()

    def handle_review_notes(self, query):
        ids = []
        raw = get_first(query, "measurement_ids", "")
        if raw:
            ids = [item for item in raw.split(",") if item]
        response_json(self, {"notes": selected_notes(ids), "statuses": STATUSES})

    def handle_save_review_note(self):
        length = parse_int(self.headers.get("Content-Length"), 0)
        body = self.rfile.read(length).decode("utf-8") if length else "{}"
        payload = json.loads(body)
        note = upsert_note(
            payload.get("measurement_id"),
            payload.get("status"),
            payload.get("note", ""),
            context=payload.get("context", {}),
        )
        response_json(self, {"note": note, "statuses": STATUSES})


HTML = r"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>DIB Spectrum Viewer</title>
  <style>
    :root {
      --bg: #f6f7f9;
      --panel: #ffffff;
      --line: #d8dde5;
      --text: #17202a;
      --muted: #5d6978;
      --accent: #1565c0;
      --primary: #c62828;
      --secondary: #6a1b9a;
      --range: rgba(21, 101, 192, 0.12);
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--text);
      background: var(--bg);
    }
    .app {
      display: grid;
      grid-template-columns: 300px minmax(360px, 1fr) 380px;
      height: 100vh;
      min-height: 620px;
    }
    aside, main, section {
      min-width: 0;
      min-height: 0;
    }
    aside {
      border-right: 1px solid var(--line);
      background: var(--panel);
      padding: 16px;
      overflow: auto;
    }
    main {
      display: grid;
      grid-template-rows: auto 1fr auto;
      padding: 16px;
      gap: 12px;
      min-width: 0;
    }
    section {
      border-left: 1px solid var(--line);
      background: var(--panel);
      padding: 16px;
      overflow: auto;
    }
    h1, h2 {
      margin: 0 0 12px;
      font-weight: 650;
      letter-spacing: 0;
    }
    h1 { font-size: 19px; }
    h2 { font-size: 15px; }
    label {
      display: block;
      margin: 10px 0 5px;
      color: var(--muted);
      font-size: 12px;
    }
    input, select, button, textarea {
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fff;
      color: var(--text);
      font-size: 13px;
      padding: 0 10px;
    }
    input, select, button { height: 34px; }
    textarea {
      min-height: 84px;
      padding: 8px 10px;
      resize: vertical;
      line-height: 1.4;
      font-family: inherit;
    }
    button {
      margin-top: 10px;
      background: var(--accent);
      border-color: var(--accent);
      color: #fff;
      cursor: pointer;
      font-weight: 600;
    }
    button.secondary {
      background: #fff;
      color: var(--accent);
    }
    .result-list, .measurement-list {
      display: grid;
      gap: 8px;
      margin-top: 12px;
    }
    .item {
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 10px;
      cursor: pointer;
      background: #fff;
    }
    .item.active {
      border-color: var(--accent);
      box-shadow: inset 3px 0 0 var(--accent);
    }
    .item-title {
      font-size: 13px;
      font-weight: 650;
      overflow-wrap: anywhere;
    }
    .item-meta, .status, .path {
      color: var(--muted);
      font-size: 12px;
      line-height: 1.45;
      overflow-wrap: anywhere;
    }
    .toolbar {
      display: grid;
      grid-template-columns: 1.3fr 0.8fr 0.8fr 0.8fr;
      gap: 10px;
      align-items: end;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 12px;
    }
    .plot-wrap {
      position: relative;
      min-height: 320px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #fff;
      overflow: hidden;
    }
    canvas {
      display: block;
      width: 100%;
      height: 100%;
    }
    .plot-status {
      position: absolute;
      left: 14px;
      top: 12px;
      color: var(--muted);
      font-size: 12px;
      background: rgba(255,255,255,0.88);
      padding: 4px 6px;
      border-radius: 4px;
      pointer-events: none;
    }
    .details {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 12px;
      display: grid;
      gap: 4px;
    }
    .grid {
      display: grid;
      grid-template-columns: 100px 1fr;
      gap: 4px 10px;
      font-size: 12px;
    }
    .key { color: var(--muted); }
    .value { overflow-wrap: anywhere; }
    .pill {
      display: inline-block;
      border-radius: 999px;
      padding: 2px 8px;
      font-size: 11px;
      background: #edf2f7;
      color: #1f2937;
      margin-right: 4px;
    }
    .pill.primary {
      background: #ffebee;
      color: var(--primary);
      font-weight: 700;
    }
    .pill.auto { background: #e8f5e9; color: #1b5e20; }
    .pill.manual { background: #fff3e0; color: #8a4b08; }
    .pill.norm { background: #e3f2fd; color: #0d47a1; }
    .pill.multi { background: #f3e5f5; color: #4a148c; }
    .pill.warn { background: #fff8e1; color: #8a5a00; font-weight: 700; }
    .pill.comment { background: #eceff1; color: #263238; }
    .pill.note { background: #e0f2f1; color: #004d40; font-weight: 700; }
    .item.suspect {
      border-color: #d7a520;
      box-shadow: inset 3px 0 0 #d7a520;
    }
    .review-editor {
      margin-top: 16px;
      border-top: 1px solid var(--line);
      padding-top: 14px;
    }
    .side-block {
      margin-top: 16px;
      border-top: 1px solid var(--line);
      padding-top: 14px;
    }
    @media (max-width: 1050px) {
      .app {
        grid-template-columns: 260px 1fr;
        grid-template-rows: 62vh 38vh;
      }
      main { grid-column: 2; grid-row: 1 / 3; }
      section { grid-column: 1; border-left: 0; border-top: 1px solid var(--line); }
    }
  </style>
</head>
<body>
  <div class="app">
    <aside>
      <h1>DIB Spectrum Viewer</h1>
      <label for="searchBox">Object search</label>
      <input id="searchBox" value="HD147889" autocomplete="off">
      <button id="searchButton">Search</button>
      <div id="objectResults" class="result-list"></div>
      <div class="side-block">
        <h2>Compare</h2>
        <label for="compareObjectIds">Object IDs</label>
        <textarea id="compareObjectIds">59,66,9</textarea>
        <label for="compareObjectId">Single objectID</label>
        <input id="compareObjectId" value="59">
        <label for="compareWavelength">Center wavelength</label>
        <input id="compareWavelength" value="13175.9">
        <label for="compareWindow">Window</label>
        <input id="compareWindow" value="25">
        <button id="compareButton">Compare Spectra</button>
        <button id="compareCombinesButton" class="secondary">Compare Combines</button>
        <button id="singleViewButton" class="secondary">Single View</button>
      </div>
    </aside>

    <main>
      <div class="toolbar">
        <div>
          <label for="combineSelect">Combine</label>
          <select id="combineSelect"></select>
        </div>
        <div>
          <label for="orderSelect">Order</label>
          <select id="orderSelect"></select>
        </div>
        <div>
          <label for="waveMin">Wave min</label>
          <input id="waveMin" placeholder="auto">
        </div>
        <div>
          <label for="waveMax">Wave max</label>
          <input id="waveMax" placeholder="auto">
        </div>
      </div>

      <div class="plot-wrap" id="plotWrap">
        <canvas id="plotCanvas"></canvas>
        <div class="plot-status" id="plotStatus">No spectrum loaded</div>
      </div>

      <div class="details">
        <h2>Selection</h2>
        <div class="grid" id="selectionDetails"></div>
      </div>
    </main>

    <section>
      <h2>Measurements</h2>
      <div class="status" id="measurementStatus">No measurements loaded</div>
      <div id="measurements" class="measurement-list"></div>
      <div class="review-editor">
        <h2>Review Note</h2>
        <div class="status" id="reviewStatus">Select a measurement</div>
        <label for="reviewStatusSelect">Status</label>
        <select id="reviewStatusSelect"></select>
        <label for="reviewText">Note</label>
        <textarea id="reviewText" placeholder="Local note only. This does not write to MySQL."></textarea>
        <button id="saveReviewButton">Save Note</button>
      </div>
    </section>
  </div>

  <script>
    const state = {
      objects: [],
      object: null,
      combines: [],
      orders: [],
      spectrumPayload: null,
      selectedMeasurementId: null,
      viewRange: {min: null, max: null, mode: 'auto'},
      reviewNotes: {},
      reviewStatuses: ['ok', 'check', 'bad_continuum', 'telluric', 'remeasure', 'ignore'],
      mode: 'single',
      comparePayload: null
    };

    const els = {
      searchBox: document.getElementById('searchBox'),
      searchButton: document.getElementById('searchButton'),
      objectResults: document.getElementById('objectResults'),
      combineSelect: document.getElementById('combineSelect'),
      orderSelect: document.getElementById('orderSelect'),
      waveMin: document.getElementById('waveMin'),
      waveMax: document.getElementById('waveMax'),
      canvas: document.getElementById('plotCanvas'),
      plotWrap: document.getElementById('plotWrap'),
      plotStatus: document.getElementById('plotStatus'),
      selectionDetails: document.getElementById('selectionDetails'),
      measurementStatus: document.getElementById('measurementStatus'),
      measurements: document.getElementById('measurements'),
      reviewStatus: document.getElementById('reviewStatus'),
      reviewStatusSelect: document.getElementById('reviewStatusSelect'),
      reviewText: document.getElementById('reviewText'),
      saveReviewButton: document.getElementById('saveReviewButton'),
      compareObjectIds: document.getElementById('compareObjectIds'),
      compareObjectId: document.getElementById('compareObjectId'),
      compareWavelength: document.getElementById('compareWavelength'),
      compareWindow: document.getElementById('compareWindow'),
      compareButton: document.getElementById('compareButton'),
      compareCombinesButton: document.getElementById('compareCombinesButton'),
      singleViewButton: document.getElementById('singleViewButton')
    };

    async function fetchJson(url) {
      const response = await fetch(url);
      const payload = await response.json();
      if (!response.ok || payload.error) {
        throw new Error(payload.error || response.statusText);
      }
      return payload;
    }

    function escapeHtml(value) {
      return String(value ?? '').replace(/[&<>"']/g, ch => ({
        '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#039;'
      }[ch]));
    }

    function renderKeyValue(container, rows) {
      container.innerHTML = rows.map(([key, value]) => `
        <div class="key">${escapeHtml(key)}</div>
        <div class="value">${escapeHtml(value ?? '')}</div>
      `).join('');
    }

    function numberOrNull(value) {
      const num = Number(value);
      return Number.isFinite(num) ? num : null;
    }

    function setWaveInputs(minValue, maxValue) {
      els.waveMin.value = minValue === null || minValue === undefined ? '' : Number(minValue).toFixed(2);
      els.waveMax.value = maxValue === null || maxValue === undefined ? '' : Number(maxValue).toFixed(2);
    }

    function captureWaveRange(mode) {
      const minValue = numberOrNull(els.waveMin.value.trim());
      const maxValue = numberOrNull(els.waveMax.value.trim());
      if (minValue !== null && maxValue !== null && maxValue > minValue) {
        state.viewRange = {min: minValue, max: maxValue, mode: mode || 'manual'};
      } else {
        state.viewRange = {min: null, max: null, mode: 'auto'};
      }
    }

    function setMeasurementRange(row) {
      const start = numberOrNull(row.integration_start);
      const end = numberOrNull(row.integration_end);
      const center = numberOrNull(row.centerlam_air) ?? numberOrNull(row.wavelength_air);
      if (start !== null && end !== null && end > start) {
        const width = Math.max(end - start, 12);
        const pad = Math.max(width * 0.75, 6);
        state.viewRange = {min: start - pad, max: end + pad, mode: 'measurement'};
      } else if (center !== null) {
        state.viewRange = {min: center - 18, max: center + 18, mode: 'measurement'};
      } else {
        state.viewRange = {min: null, max: null, mode: 'auto'};
      }
      setWaveInputs(state.viewRange.min, state.viewRange.max);
    }

    function selectBestOrderForRange(orders, currentOrder) {
      if (!orders.length) return null;
      const current = orders.find(row => String(row.echelle_order) === String(currentOrder));
      if (current) return current;
      const minValue = state.viewRange.min;
      const maxValue = state.viewRange.max;
      if (minValue !== null && maxValue !== null) {
        const center = (minValue + maxValue) / 2;
        const covering = orders.find(row => Number(row.lambdamin) <= center && Number(row.lambdamax) >= center);
        if (covering) return covering;
      }
      return orders.find(row => Number(row.echelle_order) === 42) || orders[0];
    }

    function selectedMeasurement() {
      const rows = (state.spectrumPayload && state.spectrumPayload.measurements) || [];
      return rows.find(row => row.measurement_id === state.selectedMeasurementId) || null;
    }

    function hasText(value) {
      return value !== null && value !== undefined && String(value).trim() !== '' && String(value).trim() !== 'INDEF';
    }

    function measurementFlags(row) {
      const flags = [];
      if (Number(row.primaryflag) === 1) flags.push({label: 'primary', type: 'primary'});
      flags.push(Number(row.automeasurementflag) === 1 ? {label: 'auto', type: 'auto'} : {label: 'manual', type: 'manual'});
      if (Number(row.autonormalizeflag) === 1) flags.push({label: 'norm-auto', type: 'norm'});
      if (Number(row.multiorderflag) === 1) flags.push({label: 'multi-order', type: 'multi'});
      if (numberOrNull(row.ew) === 0) flags.push({label: 'EW=0', type: 'warn'});
      if (numberOrNull(row.fwhm) === 0) flags.push({label: 'FWHM=0', type: 'warn'});
      if (numberOrNull(row.snr) !== null && numberOrNull(row.snr) < 50) flags.push({label: 'low S/N', type: 'warn'});
      if (hasText(row.comment)) flags.push({label: 'comment', type: 'comment'});
      return flags;
    }

    function isSuspectMeasurement(row) {
      return measurementFlags(row).some(flag => flag.type === 'warn');
    }

    function renderPills(flags) {
      return flags.map(flag => `<span class="pill ${flag.type}">${escapeHtml(flag.label)}</span>`).join('');
    }

    function selectedReviewNote() {
      if (!state.selectedMeasurementId) return null;
      return state.reviewNotes[state.selectedMeasurementId] || null;
    }

    function renderReviewEditor() {
      els.reviewStatusSelect.innerHTML = state.reviewStatuses.map(status => `
        <option value="${escapeHtml(status)}">${escapeHtml(status)}</option>
      `).join('');
      const measurement = selectedMeasurement();
      if (!measurement) {
        els.reviewStatus.textContent = 'Select a measurement';
        els.reviewStatusSelect.value = 'check';
        els.reviewText.value = '';
        els.reviewStatusSelect.disabled = true;
        els.reviewText.disabled = true;
        els.saveReviewButton.disabled = true;
        return;
      }

      const note = selectedReviewNote();
      els.reviewStatus.textContent = note ? `Saved: ${note.updated_at}` : `No note for ${measurement.measurement_id}`;
      els.reviewStatusSelect.disabled = false;
      els.reviewText.disabled = false;
      els.saveReviewButton.disabled = false;
      els.reviewStatusSelect.value = note ? note.status : 'check';
      els.reviewText.value = note ? note.note : '';
    }

    async function loadReviewNotes() {
      const rows = (state.spectrumPayload && state.spectrumPayload.measurements) || [];
      const ids = rows.map(row => row.measurement_id).join(',');
      if (!ids) {
        state.reviewNotes = {};
        renderReviewEditor();
        return;
      }
      const payload = await fetchJson(`/api/review-notes?measurement_ids=${encodeURIComponent(ids)}`);
      state.reviewNotes = payload.notes || {};
      state.reviewStatuses = payload.statuses || state.reviewStatuses;
      renderMeasurements();
      renderReviewEditor();
    }

    async function saveReviewNote() {
      const measurement = selectedMeasurement();
      if (!measurement) return;
      const payload = {
        measurement_id: measurement.measurement_id,
        status: els.reviewStatusSelect.value,
        note: els.reviewText.value,
        context: {
          object_id: state.object ? state.object.object_id : null,
          objectname: state.object ? state.object.objectname : null,
          combine_id: els.combineSelect.value || null,
          order: els.orderSelect.value || null,
          dib_id: measurement.dib_id,
          wavelength_air: measurement.wavelength_air
        }
      };
      const response = await fetch('/api/review-note', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify(payload)
      });
      const saved = await response.json();
      if (!response.ok || saved.error) {
        throw new Error(saved.error || response.statusText);
      }
      state.reviewNotes[saved.note.measurement_id] = saved.note;
      state.reviewStatuses = saved.statuses || state.reviewStatuses;
      renderSelection();
      renderMeasurements();
      renderReviewEditor();
    }

    async function searchObjects() {
      const q = encodeURIComponent(els.searchBox.value.trim());
      els.objectResults.innerHTML = '<div class="status">Searching...</div>';
      const rows = await fetchJson(`/api/objects?q=${q}&limit=25`);
      state.objects = rows;
      els.objectResults.innerHTML = rows.map(row => `
        <div class="item" data-object-id="${row.object_id}">
          <div class="item-title">${escapeHtml(row.objectname)} <span class="pill">ID ${row.object_id}</span></div>
          <div class="item-meta">${escapeHtml(row.aliases || '')}</div>
          <div class="item-meta">${escapeHtml(row.sptype || '')}  E(B-V): ${escapeHtml(row.e_bv)}</div>
        </div>
      `).join('') || '<div class="status">No objects found</div>';
    }

    async function selectObject(objectId) {
      [...els.objectResults.querySelectorAll('.item')].forEach(node => {
        node.classList.toggle('active', node.dataset.objectId === String(objectId));
      });
      const payload = await fetchJson(`/api/object?object_id=${encodeURIComponent(objectId)}`);
      state.object = payload.object;
      state.combines = payload.combines;
      state.viewRange = {min: null, max: null, mode: 'auto'};
      setWaveInputs(null, null);
      els.combineSelect.innerHTML = state.combines.map(row => `
        <option value="${escapeHtml(row.combine_id)}">${escapeHtml(row.combine_id)}</option>
      `).join('');
      renderSelection();
      await loadOrdersAndSpectrum();
    }

    async function loadOrdersAndSpectrum() {
      const combineId = els.combineSelect.value;
      if (!combineId) return;
      const previousOrder = els.orderSelect.value;
      const allOrders = await fetchJson(`/api/orders?combine_id=${encodeURIComponent(combineId)}`);
      state.orders = allOrders;
      els.orderSelect.innerHTML = allOrders.map(row => `
        <option value="${row.echelle_order}">${row.echelle_order}</option>
      `).join('');
      const preferred = selectBestOrderForRange(allOrders, previousOrder);
      if (preferred) {
        els.orderSelect.value = String(preferred.echelle_order);
        await loadSpectrum();
      }
    }

    async function loadSpectrum() {
      state.mode = 'single';
      const combineId = els.combineSelect.value;
      const order = els.orderSelect.value;
      if (!combineId || !order) return;

      const params = new URLSearchParams({
        combine_id: combineId,
        order: order,
        max_points: '3000'
      });
      if (state.viewRange.min !== null && state.viewRange.max !== null) {
        params.set('wave_min', state.viewRange.min);
        params.set('wave_max', state.viewRange.max);
      }

      els.plotStatus.textContent = 'Loading spectrum...';
      state.spectrumPayload = await fetchJson(`/api/spectrum?${params.toString()}`);
      renderSelection();
      renderMeasurements();
      await loadReviewNotes();
      drawPlot();
    }

    async function loadComparison() {
      const objectIds = els.compareObjectIds.value.trim();
      const wavelength = els.compareWavelength.value.trim();
      const windowValue = els.compareWindow.value.trim() || '25';
      if (!objectIds || !wavelength) return;
      const params = new URLSearchParams({
        object_ids: objectIds,
        wavelength: wavelength,
        window: windowValue,
        max_points: '900'
      });
      els.plotStatus.textContent = 'Loading comparison...';
      state.mode = 'compare';
      state.comparePayload = await fetchJson(`/api/compare?${params.toString()}`);
      renderComparisonSelection();
      drawComparisonPlot();
    }

    async function loadCombineComparison() {
      const objectId = els.compareObjectId.value.trim();
      const wavelength = els.compareWavelength.value.trim();
      const windowValue = els.compareWindow.value.trim() || '25';
      if (!objectId || !wavelength) return;
      const params = new URLSearchParams({
        object_id: objectId,
        wavelength: wavelength,
        window: windowValue,
        max_points: '900'
      });
      els.plotStatus.textContent = 'Loading combine comparison...';
      state.mode = 'compare';
      state.comparePayload = await fetchJson(`/api/compare-combines?${params.toString()}`);
      renderComparisonSelection();
      drawComparisonPlot();
    }

    function renderComparisonSelection() {
      const payload = state.comparePayload || {};
      const rows = (payload.rows || []);
      renderKeyValue(els.selectionDetails, [
        ['Mode', payload.comparison_type === 'combines' ? 'combine comparison' : 'object comparison'],
        ['Object', payload.object ? `${payload.object.objectname} / ID ${payload.object.object_id}` : ''],
        ['Objects', (payload.object_ids || []).join(', ')],
        ['Center', payload.wavelength || ''],
        ['Range', payload.wave_min ? `${Number(payload.wave_min).toFixed(2)} - ${Number(payload.wave_max).toFixed(2)}` : ''],
        ['Rows', `${rows.length}`],
        ['Errors', `${rows.filter(row => row.error).length}`]
      ]);
      els.measurementStatus.textContent = 'Comparison mode';
      els.measurements.innerHTML = rows.map(row => {
        if (row.error) {
          const label = row.object ? row.object.objectname : `objectID ${row.object_id}`;
          return `<div class="item suspect"><div class="item-title">${escapeHtml(label)}</div><div class="item-meta">${escapeHtml(row.error)}</div></div>`;
        }
        const near = (row.measurements || []).filter(item => (
          Number(item.wavelength_air) >= Number(payload.wave_min) &&
          Number(item.wavelength_air) <= Number(payload.wave_max)
        ));
        return `
          <div class="item">
            <div class="item-title">${escapeHtml(row.object.objectname)} <span class="pill">ID ${row.object.object_id}</span></div>
            <div class="item-meta">${escapeHtml(row.combine.combine_id)} / order ${row.order.echelle_order}</div>
            <div class="item-meta">telluric ${row.combine.telluricflag}, combine ${row.combine.combineflag}, number ${row.combine.combine_number}</div>
            <div class="item-meta">${near.length} measurements in range</div>
          </div>
        `;
      }).join('');
      renderReviewEditor();
    }

    function renderSelection() {
      const obj = state.object || {};
      const payload = state.spectrumPayload || {};
      const file = payload.spectrum_file || {};
      const measurement = selectedMeasurement();
      const rows = [
        ['Object', obj.objectname ? `${obj.objectname} / ID ${obj.object_id}` : ''],
        ['Aliases', obj.aliases || ''],
        ['Type', obj.sptype || ''],
        ['Coordinates', obj.ra && obj.decli ? `${obj.ra} ${obj.decli}` : ''],
        ['Combine', els.combineSelect.value || ''],
        ['Order', els.orderSelect.value || ''],
        ['View', state.viewRange.mode === 'auto' ? 'auto' : `${state.viewRange.min.toFixed(2)} - ${state.viewRange.max.toFixed(2)} (${state.viewRange.mode})`],
        ['Range', file.lambdamin ? `${Number(file.lambdamin).toFixed(2)} - ${Number(file.lambdamax).toFixed(2)}` : ''],
        ['FITS', file.combinefilepath || '']
      ];
      if (measurement) {
        const note = selectedReviewNote();
        rows.push(
          ['Measurement', measurement.measurement_id],
          ['DIB', `${measurement.dib_id} / ${Number(measurement.wavelength_air).toFixed(2)} A / ${measurement.category || ''}`],
          ['Flags', measurementFlags(measurement).map(flag => flag.label).join(', ')],
          ['Review', note ? `${note.status}: ${note.note}` : ''],
          ['EW', `${measurement.ew} +/- ${measurement.ew_err}`],
          ['Center', `${measurement.centerlam_air} A / v=${measurement.helio_velocity}`],
          ['FWHM/SNR', `${measurement.fwhm} / ${measurement.snr}`],
          ['Integration', `${measurement.integration_start} - ${measurement.integration_end}`],
          ['Comment', measurement.comment || '']
        );
      }
      renderKeyValue(els.selectionDetails, rows);
    }

    function renderMeasurements() {
      const rows = (state.spectrumPayload && state.spectrumPayload.measurements) || [];
      const suspectCount = rows.filter(isSuspectMeasurement).length;
      els.measurementStatus.textContent = `${rows.length} measurements / ${suspectCount} flagged`;
      els.measurements.innerHTML = rows.map(row => {
        const flags = measurementFlags(row);
        const classes = [
          'item',
          row.measurement_id === state.selectedMeasurementId ? 'active' : '',
          isSuspectMeasurement(row) ? 'suspect' : ''
        ].filter(Boolean).join(' ');
        return `
        <div class="${classes}" data-measurement-id="${escapeHtml(row.measurement_id)}">
          <div class="item-title">
            DIB ${Number(row.wavelength_air).toFixed(1)}
            <span class="pill">ID ${row.dib_id}</span>
            ${renderPills(flags)}
            ${state.reviewNotes[row.measurement_id] ? `<span class="pill note">${escapeHtml(state.reviewNotes[row.measurement_id].status)}</span>` : ''}
          </div>
          <div class="item-meta">EW ${row.ew} +/- ${row.ew_err}, FWHM ${row.fwhm}</div>
          <div class="item-meta">center ${row.centerlam_air}, S/N ${row.snr}, range ${row.integration_start} - ${row.integration_end}</div>
        </div>
      `;
      }).join('') || '<div class="status">No measurements for this order</div>';
    }

    function getPlotBounds(x, y) {
      const xmin = Math.min(...x);
      const xmax = Math.max(...x);
      const finiteY = y.filter(Number.isFinite);
      const sortedY = [...finiteY].sort((a, b) => a - b);
      const q = p => sortedY[Math.max(0, Math.min(sortedY.length - 1, Math.floor((sortedY.length - 1) * p)))];
      let ymin = q(0.01);
      let ymax = q(0.99);
      if (ymin === ymax) {
        ymin -= 0.05;
        ymax += 0.05;
      }
      const pad = (ymax - ymin) * 0.08;
      return {xmin, xmax, ymin: ymin - pad, ymax: ymax + pad};
    }

    function drawPlot() {
      if (state.mode === 'compare') {
        drawComparisonPlot();
        return;
      }
      const payload = state.spectrumPayload;
      const canvas = els.canvas;
      const rect = els.plotWrap.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      canvas.width = Math.max(320, Math.floor(rect.width * dpr));
      canvas.height = Math.max(260, Math.floor(rect.height * dpr));
      const ctx = canvas.getContext('2d');
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      const width = canvas.width / dpr;
      const height = canvas.height / dpr;
      ctx.clearRect(0, 0, width, height);

      if (!payload || !payload.spectrum || payload.spectrum.wavelength.length < 2) {
        els.plotStatus.textContent = 'No spectrum loaded';
        return;
      }

      const x = payload.spectrum.wavelength;
      const y = payload.spectrum.flux;
      const bounds = getPlotBounds(x, y);
      const margin = {left: 58, right: 18, top: 20, bottom: 46};
      const pw = width - margin.left - margin.right;
      const ph = height - margin.top - margin.bottom;
      const sx = v => margin.left + (v - bounds.xmin) / (bounds.xmax - bounds.xmin) * pw;
      const sy = v => margin.top + (bounds.ymax - v) / (bounds.ymax - bounds.ymin) * ph;

      ctx.fillStyle = '#fff';
      ctx.fillRect(0, 0, width, height);
      ctx.strokeStyle = '#d8dde5';
      ctx.lineWidth = 1;
      ctx.strokeRect(margin.left, margin.top, pw, ph);

      ctx.fillStyle = '#5d6978';
      ctx.font = '12px -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif';
      ctx.textAlign = 'center';
      for (let i = 0; i <= 4; i++) {
        const xv = bounds.xmin + (bounds.xmax - bounds.xmin) * i / 4;
        const px = sx(xv);
        ctx.strokeStyle = '#eef1f5';
        ctx.beginPath();
        ctx.moveTo(px, margin.top);
        ctx.lineTo(px, margin.top + ph);
        ctx.stroke();
        ctx.fillText(xv.toFixed(1), px, height - 18);
      }
      ctx.textAlign = 'right';
      for (let i = 0; i <= 4; i++) {
        const yv = bounds.ymin + (bounds.ymax - bounds.ymin) * i / 4;
        const py = sy(yv);
        ctx.strokeStyle = '#eef1f5';
        ctx.beginPath();
        ctx.moveTo(margin.left, py);
        ctx.lineTo(margin.left + pw, py);
        ctx.stroke();
        ctx.fillText(yv.toFixed(2), margin.left - 8, py + 4);
      }

      for (const row of payload.measurements || []) {
        if (row.integration_start && row.integration_end) {
          const x0 = Math.max(margin.left, sx(row.integration_start));
          const x1 = Math.min(margin.left + pw, sx(row.integration_end));
          if (x1 > margin.left && x0 < margin.left + pw) {
            ctx.fillStyle = row.primaryflag ? 'rgba(198, 40, 40, 0.12)' : 'rgba(21, 101, 192, 0.10)';
            ctx.fillRect(x0, margin.top, x1 - x0, ph);
          }
        }
      }

      ctx.strokeStyle = '#17202a';
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      for (let i = 0; i < x.length; i++) {
        const px = sx(x[i]);
        const py = sy(y[i]);
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.stroke();

      for (const row of payload.measurements || []) {
        const color = row.primaryflag ? '#c62828' : '#6a1b9a';
        if (row.wavelength_air >= bounds.xmin && row.wavelength_air <= bounds.xmax) {
          ctx.strokeStyle = color;
          ctx.setLineDash([5, 4]);
          ctx.beginPath();
          ctx.moveTo(sx(row.wavelength_air), margin.top);
          ctx.lineTo(sx(row.wavelength_air), margin.top + ph);
          ctx.stroke();
          ctx.setLineDash([]);
        }
        if (row.centerlam_air >= bounds.xmin && row.centerlam_air <= bounds.xmax) {
          ctx.strokeStyle = color;
          ctx.lineWidth = 2;
          ctx.beginPath();
          ctx.moveTo(sx(row.centerlam_air), margin.top);
          ctx.lineTo(sx(row.centerlam_air), margin.top + ph);
          ctx.stroke();
          ctx.lineWidth = 1;
        }
      }

      els.plotStatus.textContent = `${payload.combine_id}  order ${payload.order}  ${payload.spectrum.points}/${payload.spectrum.full_points} points`;
    }

    function drawComparisonPlot() {
      const payload = state.comparePayload;
      const canvas = els.canvas;
      const rect = els.plotWrap.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      canvas.width = Math.max(320, Math.floor(rect.width * dpr));
      canvas.height = Math.max(260, Math.floor(rect.height * dpr));
      const ctx = canvas.getContext('2d');
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      const width = canvas.width / dpr;
      const height = canvas.height / dpr;
      ctx.clearRect(0, 0, width, height);
      ctx.fillStyle = '#fff';
      ctx.fillRect(0, 0, width, height);

      if (!payload || !payload.rows || !payload.rows.length) {
        els.plotStatus.textContent = 'No comparison loaded';
        return;
      }

      const rows = payload.rows.filter(row => row.spectrum && row.spectrum.wavelength.length > 1);
      if (!rows.length) {
        els.plotStatus.textContent = 'No comparable spectra';
        return;
      }

      const margin = {left: 76, right: 18, top: 24, bottom: 36};
      const pw = width - margin.left - margin.right;
      const ph = height - margin.top - margin.bottom;
      const xmin = Number(payload.wave_min);
      const xmax = Number(payload.wave_max);
      const rowHeight = ph / rows.length;
      const sx = v => margin.left + (v - xmin) / (xmax - xmin) * pw;

      ctx.strokeStyle = '#d8dde5';
      ctx.lineWidth = 1;
      ctx.strokeRect(margin.left, margin.top, pw, ph);

      ctx.fillStyle = '#5d6978';
      ctx.font = '12px -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif';
      ctx.textAlign = 'center';
      for (let i = 0; i <= 4; i++) {
        const xv = xmin + (xmax - xmin) * i / 4;
        const px = sx(xv);
        ctx.strokeStyle = '#eef1f5';
        ctx.beginPath();
        ctx.moveTo(px, margin.top);
        ctx.lineTo(px, margin.top + ph);
        ctx.stroke();
        ctx.fillStyle = '#5d6978';
        ctx.fillText(xv.toFixed(1), px, height - 14);
      }

      const centerX = sx(Number(payload.wavelength));
      ctx.strokeStyle = '#c62828';
      ctx.setLineDash([5, 4]);
      ctx.beginPath();
      ctx.moveTo(centerX, margin.top);
      ctx.lineTo(centerX, margin.top + ph);
      ctx.stroke();
      ctx.setLineDash([]);

      rows.forEach((row, index) => {
        const x = row.spectrum.wavelength;
        const y = row.spectrum.flux;
        const bounds = getPlotBounds(x, y);
        const top = margin.top + rowHeight * index + 8;
        const bottom = margin.top + rowHeight * (index + 1) - 8;
        const rh = Math.max(18, bottom - top);
        const sy = v => top + (bounds.ymax - v) / (bounds.ymax - bounds.ymin) * rh;

        ctx.strokeStyle = '#e6eaf0';
        ctx.beginPath();
        ctx.moveTo(margin.left, margin.top + rowHeight * index);
        ctx.lineTo(margin.left + pw, margin.top + rowHeight * index);
        ctx.stroke();

        for (const measurement of row.measurements || []) {
          if (measurement.integration_start && measurement.integration_end) {
            const x0 = sx(Number(measurement.integration_start));
            const x1 = sx(Number(measurement.integration_end));
            if (x1 > margin.left && x0 < margin.left + pw) {
              ctx.fillStyle = Number(measurement.primaryflag) === 1 ? 'rgba(198, 40, 40, 0.12)' : 'rgba(21, 101, 192, 0.08)';
              ctx.fillRect(Math.max(margin.left, x0), top, Math.min(margin.left + pw, x1) - Math.max(margin.left, x0), rh);
            }
          }
        }

        ctx.strokeStyle = '#17202a';
        ctx.lineWidth = 1.1;
        ctx.beginPath();
        x.forEach((xv, pointIndex) => {
          const px = sx(Number(xv));
          const py = sy(Number(y[pointIndex]));
          if (pointIndex === 0) ctx.moveTo(px, py);
          else ctx.lineTo(px, py);
        });
        ctx.stroke();

        ctx.fillStyle = '#17202a';
        ctx.textAlign = 'right';
        const label = payload.comparison_type === 'combines'
          ? `${row.combine.combine_id} / m${row.order.echelle_order}`
          : `${row.object.objectname} / o${row.object.object_id} / m${row.order.echelle_order}`;
        ctx.fillText(label, margin.left - 8, top + 14);
      });

      const label = payload.comparison_type === 'combines' ? 'combine comparison' : 'object comparison';
      els.plotStatus.textContent = `${label}  ${rows.length}/${payload.rows.length} spectra  ${payload.wave_min.toFixed(1)}-${payload.wave_max.toFixed(1)}`;
    }

    els.searchButton.addEventListener('click', () => searchObjects().catch(err => alert(err.message)));
    els.searchBox.addEventListener('keydown', event => {
      if (event.key === 'Enter') searchObjects().catch(err => alert(err.message));
    });
    els.objectResults.addEventListener('click', event => {
      const item = event.target.closest('.item[data-object-id]');
      if (item) selectObject(item.dataset.objectId).catch(err => alert(err.message));
    });
    els.combineSelect.addEventListener('change', () => {
      captureWaveRange(state.viewRange.mode === 'auto' ? 'auto' : 'manual');
      loadOrdersAndSpectrum().catch(err => alert(err.message));
    });
    els.orderSelect.addEventListener('change', () => {
      captureWaveRange(state.viewRange.mode === 'auto' ? 'auto' : 'manual');
      loadSpectrum().catch(err => alert(err.message));
    });
    els.waveMin.addEventListener('keydown', event => {
      if (event.key === 'Enter') {
        captureWaveRange('manual');
        loadSpectrum().catch(err => alert(err.message));
      }
    });
    els.waveMax.addEventListener('keydown', event => {
      if (event.key === 'Enter') {
        captureWaveRange('manual');
        loadSpectrum().catch(err => alert(err.message));
      }
    });
    els.measurements.addEventListener('click', event => {
      const item = event.target.closest('.item[data-measurement-id]');
      if (!item) return;
      state.selectedMeasurementId = item.dataset.measurementId;
      [...els.measurements.querySelectorAll('.item')].forEach(node => {
        node.classList.toggle('active', node.dataset.measurementId === state.selectedMeasurementId);
      });
      const row = (state.spectrumPayload.measurements || []).find(m => m.measurement_id === state.selectedMeasurementId);
      if (row) {
        renderSelection();
        setMeasurementRange(row);
        renderReviewEditor();
        loadSpectrum().catch(err => alert(err.message));
      }
    });
    els.saveReviewButton.addEventListener('click', () => {
      saveReviewNote().catch(err => alert(err.message));
    });
    els.compareButton.addEventListener('click', () => {
      loadComparison().catch(err => alert(err.message));
    });
    els.compareCombinesButton.addEventListener('click', () => {
      loadCombineComparison().catch(err => alert(err.message));
    });
    els.singleViewButton.addEventListener('click', () => {
      state.mode = 'single';
      drawPlot();
      renderSelection();
      renderMeasurements();
      renderReviewEditor();
    });
    window.addEventListener('resize', () => drawPlot());

    searchObjects()
      .then(() => {
        const first = state.objects[0];
        if (first) return selectObject(first.object_id);
      })
      .catch(err => alert(err.message));
  </script>
</body>
</html>
"""


def build_parser():
    parser = argparse.ArgumentParser(description="Start the read-only DIB spectrum viewer.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    server = ThreadedHTTPServer((args.host, args.port), ViewerHandler)
    print("DIB spectrum viewer: http://{}:{}/".format(args.host, args.port))
    print("mode: read-only")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nshutting down")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
