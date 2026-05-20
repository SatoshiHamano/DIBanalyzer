# DIBproject

This repository contains scripts and data products for DIB analysis workflows.

Current organization notes are tracked in `SCRIPT_ORGANIZATION.md`.
Database notes are tracked in `DB_SCHEMA_NOTES.md`.

Important archive directories:

- `IRAF-ver-backup/`: historical IRAF-dependent scripts.
- `figureFactory_weakDIB/`: scripts and outputs used for a published weak-DIB
  paper.

The root Python scripts are currently the active workflow, but many of them use
same-directory imports.  Move them only after import paths and smoke tests are
prepared.
