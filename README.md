# DIBproject

This repository contains scripts and data products for DIB analysis workflows.

Current organization notes are tracked in `SCRIPT_ORGANIZATION.md`.
Database notes are tracked in `DB_SCHEMA_NOTES.md`.
Workflow notes are tracked in `WORKFLOW.md`.
DB-registration workflow issues are tracked in `DB_REGISTRATION_WORKFLOW.md`.
Published-data release planning is tracked in `DATA_RELEASE_PLAN.md`.

Important archive directories:

- `IRAF-ver-backup/`: historical IRAF-dependent scripts.
- `figureFactory_weakDIB/`: scripts and outputs used for a published weak-DIB
  paper.

The root Python scripts are currently the active workflow, but many of them use
same-directory imports.  Move them only after import paths and smoke tests are
prepared.
