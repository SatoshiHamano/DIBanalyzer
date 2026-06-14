#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Local JSON review notes for the read-only DIB workbench.

These notes are intentionally kept outside MySQL.  They can be inspected,
versioned, and later converted into a dry-run DB update if needed.
"""

from __future__ import print_function

import json
import os
from datetime import datetime


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_NOTES_PATH = os.path.join(REPO_ROOT, "workbench_review_notes.json")

STATUSES = [
    "ok",
    "check",
    "bad_continuum",
    "telluric",
    "remeasure",
    "ignore",
]


def utc_now():
    return datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def load_notes(path=None):
    path = path or DEFAULT_NOTES_PATH
    if not os.path.exists(path):
        return {}
    with open(path, "r") as handle:
        payload = json.load(handle)
    if isinstance(payload, dict):
        return payload
    raise ValueError("Review notes file must contain a JSON object: {}".format(path))


def save_notes(notes, path=None):
    path = path or DEFAULT_NOTES_PATH
    directory = os.path.dirname(path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)
    tmp_path = path + ".tmp"
    with open(tmp_path, "w") as handle:
        json.dump(notes, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp_path, path)


def selected_notes(measurement_ids, path=None):
    notes = load_notes(path)
    return {
        measurement_id: notes[measurement_id]
        for measurement_id in measurement_ids
        if measurement_id in notes
    }


def upsert_note(measurement_id, status, note, context=None, path=None):
    if not measurement_id:
        raise ValueError("measurement_id is required")
    if status not in STATUSES:
        raise ValueError("Unknown review status: {}".format(status))

    notes = load_notes(path)
    previous = notes.get(measurement_id, {})
    now = utc_now()
    notes[measurement_id] = {
        "measurement_id": measurement_id,
        "status": status,
        "note": note or "",
        "context": context or previous.get("context", {}),
        "created_at": previous.get("created_at", now),
        "updated_at": now,
    }
    save_notes(notes, path)
    return notes[measurement_id]
