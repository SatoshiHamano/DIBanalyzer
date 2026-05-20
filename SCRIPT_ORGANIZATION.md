# Script organization notes

This repository currently keeps many executable scripts, analysis helpers,
generated products, and historical snapshots in the project root.  The notes
below document the intended role of the major groups before moving files.

## Archive directories

- `IRAF-ver-backup/`
  - Historical scripts from the IRAF-dependent version of the workflow.
  - Keep as an archive/reference snapshot.
  - Do not treat as the active implementation unless explicitly restoring old
    IRAF behavior.

- `figureFactory_weakDIB/`
  - Scripts and outputs used for a published weak-DIB paper.
  - Keep as a publication-reproducibility archive.
  - Avoid mixing these scripts into active pipeline refactors.

## Active root scripts

The Python files in the repository root appear to be the active or most recent
workflow.  They are still tightly coupled through same-directory imports, so
they should not be moved until imports are updated or compatibility wrappers are
added.

Most active workflow scripts also depend on the MySQL database documented in
`DB_SCHEMA_NOTES.md`.

Suggested future grouping:

- Core utilities
  - `Spec1Dtools.py`
  - `open_mysql_project.py`
  - `vac2air_spec.py`
  - `spectra_plotter.py`
  - `regression.py`

- DIB measurement and summaries
  - `DIBanalysis*.py`
  - `DIBcut.py`
  - `DIBmeasurement*.py`
  - `DIBprimary*.py`
  - `DIBsummary*.py`
  - `DIBvelocity_plot.py`

- Telluric and atmospheric correction
  - `telluric_*.py`
  - `check_*telluric.py`
  - `reference_advanced_telluric.py`
  - `atran_model_npz.py`
  - `winered_atmospheric_dispersion.py`

- Waveshift and wavelength correction
  - `waveshift_*.py`
  - `Waveshift_*.py`
  - `ccwaveshift.py`
  - `dopcor_vcorr.py`

- Database maintenance
  - `*_mysql*.py`
  - `open_mysql_project.py`
  - `object_*.py`
  - `observation_*.py`
  - `update_*.py`
  - `set_primaryflag*.py`
  - `unset_primaryflag.py`

- Correlation and plotting
  - `correlation_*.py`
  - `DIBcorr_*.py`
  - `*_figure.py`
  - `*_plot.py`
  - `ebv_sptype_dist*.py`
  - `snr_distribution.py`

- Shell/list generation and inspection helpers
  - `create_shell_*.py`
  - `find_keyword.py`
  - `list_*.py`
  - `read_centersearch_log.py`

## Non-code files in the root

These should eventually be moved out of the root after confirming paths used by
the scripts.

- Generated plots: `*.png`, `*.pdf`
- FITS samples and temporary products: `*.fits`
- Tables and parameters: `*.dat`, `*.txt`, `*.csv`, `*.xls`, `*.xlsx`, `*.list`
- Runtime/cache files: `__pycache__/`, `*.pyc`, `.DS_Store`, `*.log`

## Recommended migration order

1. Add `.gitignore` rules for cache and local runtime products.
2. Move only clearly generated files first, not active scripts.
3. Keep archive directories intact and document their purpose.
4. Convert active root scripts into a package or add compatibility wrappers.
5. Update imports and run smoke tests after each group is moved.
