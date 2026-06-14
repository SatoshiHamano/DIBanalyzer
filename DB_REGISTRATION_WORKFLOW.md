# DB registration workflow and open issues

この文書は、DIBproject の現行作業を「MySQL データベースへ何を登録するか」
を軸に並べ直した運用メモです。実行順、入力、更新されるテーブル、確認方法、
使いづらさを同じ粒度で見えるようにすることを目的にしています。

関連メモ:

- `DB_SCHEMA_NOTES.md`: DB テーブルと推定リレーション
- `WORKFLOW.md`: 旧来の作業順序の全体像
- `WORKBENCH_DB_PERFORMANCE.md`: 読み取り系クエリと index 検討
- `GUI_WORKBENCH_PLAN.md`: 軽量 GUI の計画

## Core policy

- DB は外部キー制約を持たないため、`objectID`, `pipelineID`,
  `telluricID`, `combineID`, `DIBID`, `measurementID` の命名規約を
  スクリプト側で守る必要がある。
- 現行スクリプトの多くは DB 書き込みとファイル生成を同時に行う。
- まずは科学処理を全面置換せず、各ステップに dry-run/status/check を追加する。
- 書き込み処理は将来的に `dibctl.py` から呼べる薄い入口へ集約する。

## Current registration chain

| Order | Operation | Main script | Main input | Tables written |
| ---: | --- | --- | --- | --- |
| 1 | Object registration | `object_add_mysql.py`, `object_add_mysql_ver2.py` | object name list | `object`, `objectdict` |
| 2 | Observation import | `observation_mysql_update_csv.py`, `observation_mysql_update.py` | `observationData/pinot*` or WODB HTML | `observation` |
| 3 | Pipeline product import | `datadownload_merlot_ver6.py` | reduction target list / merlot paths | `datareduction`, `reducedframe` |
| 4 | Waveshift measurement/correction | `waveshift_measure.py`, `waveshift_correct.py`, `Waveshift_fit_ver3.py` | pipeline products | mostly files; DB is read |
| 5 | Telluric correction | `telluric_auto.py` | target/ref `pipelineID` | `telluriccorrection`, `telluricresult` |
| 6 | Combine spectra | `combine_MySQL.py` | `pipelineID` or `telluricID` | `combinesummary`, `combinedataset`, `combinedspectrum` |
| 7 | DIB measurement | `DIBanalysis.py` | `combineID`, `DIBID` | `DIBmeasurement` |
| 8 | Primary selection | `set_primaryflag.py`, `set_primaryflag_list.py` | `measurementID` | `DIBmeasurement.primaryflag` |
| 9 | Literature EW import | `DIBmeasurement_literature.py` | tab-separated literature file | `DIBEWsummary` |
| 10 | Summary/figures/correlation | `DIBprimary.py`, `DIBsummary_*.py`, `DIBcorr_*.py` | object/DIB lists | mostly DB read and files |

## Step details

### 1. Object registration

Purpose:

- Register astronomical targets and telluric standards.
- Maintain alias lookup through `objectdict`.

Known scripts:

- `object_add_mysql.py`
- `object_add_mysql_ver2.py`

DB writes:

- `object`: canonical object row, type, coordinates, spectral type, E(B-V),
  photometry, radial velocity, comments.
- `objectdict`: aliases and priority name.

Inputs and dependencies:

- One object name per line.
- SIMBAD HTML response is parsed inside the script.
- `name#same` appears to mean "alias of the previous object".

Checks to add:

- Given an input name, show exact/alias match before inserting.
- Show existing coordinates and HD-style aliases.
- Detect possible duplicates by coordinate tolerance.

Current issues:

- `object_add_mysql.py` and `object_add_mysql_ver2.py` are both present.
- SIMBAD parsing depends on HTML line positions.
- Interactive `OBJECT`/`STANDARD` prompt blocks automation.
- No dry-run or duplicate report.
- DB writes are built with string formatting.

Improvement tasks:

- Decide the current object importer.
- Add `dibctl.py object lookup NAME`.
- Add `dibctl.py object import FILE --type OBJECT|STANDARD --dry-run`.
- Move SIMBAD parsing behind a small adapter and cache raw lookup results.

### 2. Observation import

Purpose:

- Register individual observed frames and metadata.
- Tie frame IDs to `objectID`.

Known scripts:

- `observation_mysql_update_csv.py`
- `observation_mysql_update.py`

DB writes:

- `observation`

Inputs and dependencies:

- `observationData/pinot*` CSV files for the newer local workflow.
- Hard-coded WODB HTML URL for the older workflow.
- `objectdict.registeredname` must already contain target aliases.

Checks to add:

- Count CSV rows by `DATA-TYP`, `WODBTHEM`, and `existence`.
- List unregistered targets before writing.
- Report frames already present with a different `objectID`.

Current issues:

- `observationData/` is a metadata import source, not actual FITS storage.
- Missing aliases cause import failures or warnings.
- The CSV importer silently glob-loads `observationData/pinot*`.
- The HTML importer depends on a remote page layout.

Improvement tasks:

- Add `dibctl.py observation scan-csv`.
- Add `dibctl.py observation import-csv --dry-run`.
- Create a preflight report: new frames, updated frames, unregistered objects.

### 3. Pipeline product import

Purpose:

- Download or locate reduced pipeline products.
- Register reduction-level paths and frame composition.

Known scripts:

- `datadownload_merlot_ver6.py`
- older tracked variant: `datadownload_merlot_ver5.py`

DB writes:

- `datareduction`
- `reducedframe`

Inputs and dependencies:

- Target list for pipeline data.
- SSH/SCP access to `merlot.kyoto-su.ac.jp`.
- Embedded remote root: `/media/WD_ext/PIPELINE_DATA`
- Embedded local roots:
  - `/Users/hamano/DIB_analysis/DIB_pipeline_dir`
  - `/Users/hamano/DIB_analysis/Standard_pipeline_dir`

Checks to add:

- For a `pipelineID`, verify `datareduction.path` exists.
- Verify each `reducedframe.objectframe` and `skyframe` exists in `observation`.
- Count expected order FITS products by mode.
- Compare version 5 and version 6 behavior for `badpix`/cosmic-ray logs.

Current issues:

- The script combines download, metadata parsing, and DB writes.
- SSH password flow is interactive.
- Absolute paths are embedded.
- `ver5` and `ver6` coexist; `ver6` adds cosmic-ray log handling.
- No dry-run/import-only mode.

Improvement tasks:

- Treat `datadownload_merlot_ver6.py` as the likely current importer.
- Add a read-only validator for existing `datareduction` and `reducedframe`.
- Split future command modes:
  - download only
  - inspect local pipeline directory
  - register DB rows

### 4. Waveshift measurement and correction

Purpose:

- Measure wavelength shift from telluric absorption lines.
- Produce `waveshift_measure/` and `waveshift_correct/` products used by later
  telluric correction.

Known scripts:

- `waveshift_measure.py`
- `waveshift_correct.py`
- `Waveshift_fit_ver3.py`
- archived old workflow: `legacy/waveshift_main_ver4/`

DB writes:

- Mainly reads `datareduction`; output is file products.

Inputs and dependencies:

- Pipeline product paths from `datareduction.path`.
- Telluric line list now located under
  `legacy/waveshift_main_ver4/telluric_single_list_ordered_selected3.dat`.
- `waveshift_parameter/*.py` contains historical diagnostic scripts.

Checks to add:

- For each `pipelineID`, report whether `waveshift_measure/` exists.
- Verify expected order files and per-order `.txt` measurements.
- Show wavelength-shift quality metrics before telluric correction.

Current issues:

- This is known to be a problem area.
- Multiple historical versions exist.
- Some behavior depends on exact directory naming under pipeline products.
- It is hard to see whether a spectrum has valid waveshift correction.

Improvement tasks:

- Build `dibctl.py waveshift status PIPELINE_ID`.
- Build a small comparison report using `waveshift_parameter` diagnostics.
- Make telluric correction fail early with a clear waveshift-status report.
- Decide whether line-list data should move out of `legacy/` into
  `reference_data/`.

### 5. Telluric correction

Purpose:

- Pair an object pipeline product with a telluric standard.
- Write corrected spectra and register per-order results.

Known scripts:

- `telluric_auto.py`
- related inspection: `check_advanced_telluric.py`,
  `reference_advanced_telluric.py`

DB writes:

- `telluriccorrection`
- `telluricresult`

Inputs and dependencies:

- Target `pipelineID`.
- Reference/telluric `pipelineID`.
- Existing `waveshift_measure/` products for the target.
- Pipeline path layout and order definitions.

Checks to add:

- Show all existing telluric corrections for a target.
- Compare target/ref slit and mode before writing.
- Verify every `telluricresult.telluricfilepath` exists.
- Plot before/after correction in the spectrum viewer.

Current issues:

- Interactive duplicate-run confirmation.
- Automatic, manual, advanced, and simple modes are mixed in one script.
- Existing Astropy/specutils replacement may have flux-scaling issues.
- Outputs and DB rows can diverge if a run fails partway.

Improvement tasks:

- Add `dibctl.py telluric status PIPELINE_ID`.
- Add dry-run pairing check.
- Add transaction boundary or rollback strategy for DB inserts.
- Keep the Astropy telluric validation script as a regression check.

### 6. Combine spectra

Purpose:

- Select one or more corrected or uncorrected datasets.
- Combine or register per-order spectra and create the `combineID` used by DIB
  measurement.

Known scripts:

- `combine_MySQL.py`

DB writes:

- `combinesummary`
- `combinedataset`
- `combinedspectrum`

Inputs and dependencies:

- One or more `telluricID` or `pipelineID` values.
- Existing FITS products on disk.
- Interactive frame and weight choices.

Checks to add:

- For an object, list candidate `pipelineID`, `telluricID`, and existing
  `combineID`.
- Verify all `combinedspectrum.combinefilepath` values exist.
- Confirm mode/order compatibility before writing.

Current issues:

- Interactive frame selection and weights make the result hard to reproduce.
- One script handles both single-dataset registration and multi-dataset combine.
- `telluric_spectra.pickle` is generated as a side product.
- Output roots are absolute.

Improvement tasks:

- Add `dibctl.py combine status OBJECT_ID`.
- Add non-interactive options for frame and weights.
- Save combine provenance in a machine-readable manifest next to output FITS.

### 7. DIB measurement

Purpose:

- Measure one DIB in one or more combined spectra.
- Register per-measurement result and generated plots/spectra.

Known scripts:

- `DIBanalysis.py`
- variants: `DIBanalysis_intrange.py`, `DIBanalysis_upperlimit.py`
- status: `DIBmeasurement_status.py`

DB writes:

- `DIBmeasurement`

Inputs and dependencies:

- `combineID`
- `DIBID`
- `DIBlist`
- `combinedspectrum`
- telluric spectra and masks

Checks to add:

- For `objectID`, show measured/unmeasured DIBs.
- For `combineID`, list available orders and wavelength coverage.
- Show measurements without primary flag and duplicated primary flags.

Current issues:

- Measurement can involve manual choices.
- It is easy to lose track of which DIBs are complete.
- Upper-limit and integration-range variants are separate scripts.
- Result files and DB rows can be hard to audit together.

Improvement tasks:

- Promote `DIBmeasurement_status.py` into `dibctl.py dib status OBJECT_ID`.
- Add viewer support for measurement review and primary selection candidates.
- Add a per-object "next DIBs to measure" report.

### 8. Primary selection

Purpose:

- Choose the representative measurement for each object/DIB.

Known scripts:

- `set_primaryflag.py`
- `set_primaryflag_list.py`
- inspection: `DIBprimary.py`

DB writes:

- `DIBmeasurement.primaryflag`

Inputs and dependencies:

- `measurementID`
- `DIBmeasurement` joined to `combinesummary` to find competing rows.

Checks to add:

- Show all competing measurements before update.
- Confirm exactly one primary per object/DIB after update.
- Detect missing or duplicate primary flags globally.

Current issues:

- Direct DB mutation with no dry-run.
- Batch script commits inside the loop.
- No integrated visual review before setting primary.

Improvement tasks:

- Add `dibctl.py dib primary-preview MEASUREMENT_ID`.
- Add `dibctl.py dib set-primary MEASUREMENT_ID --dry-run`.
- Move primary review into the spectrum viewer after read-only preview exists.

### 9. Literature EW import

Purpose:

- Import external published DIB equivalent-width measurements for comparison.

Known scripts:

- `DIBmeasurement_literature.py`

DB writes:

- `DIBEWsummary`

Inputs and dependencies:

- Tab-separated literature files.
- `DIBID` and `objectID` must already match DB rows.
- Current committed inputs include `literature_data/` and Hobbs DIB lists.

Checks to add:

- Validate every input `DIBID` and `objectID` before insert.
- Report duplicate rows that `INSERT IGNORE` would skip.
- Summarize by reference, DIBID, object count, and upper-limit count.

Current issues:

- Input file path and reference mapping are external to the script.
- `INSERT IGNORE` can hide malformed duplicate assumptions.
- Spectrograph/resolution mappings are hard-coded.

Improvement tasks:

- Add `dibctl.py literature validate FILE`.
- Add an explicit import manifest per literature source.
- Prefer a unique-key report over silent `INSERT IGNORE`.

### 10. Summary, figure, and correlation outputs

Purpose:

- Read DB measurements and produce paper figures, summaries, and correlation
  tables.

Known scripts:

- `DIBprimary.py`
- `DIBsummary_*.py`
- `DIBcorr_*.py`
- `correlation_*.py`

DB writes:

- Mostly none; primarily reads `DIBmeasurement`, `DIBEWsummary`,
  `combinesummary`, `object`, `DIBlist`.

Inputs and dependencies:

- Object lists are often embedded directly in scripts.
- Output files go to local analysis directories such as `Correlation/`.

Checks to add:

- Before plotting, report number of objects, measurements, missing primary
  flags, and upper limits.
- Keep generated outputs ignored unless they are publication-release artifacts.

Current issues:

- Many scripts overlap.
- Some files contain local uncommented/commented analysis edits.
- Output directories and object lists are not centralized.

Improvement tasks:

- Create a small configuration file for object/DIB selections.
- Keep publication-specific scripts archived separately from active analysis.
- Build read-only summaries in the viewer before regenerating static plots.

## Cross-cutting issues

### DB safety

- No foreign keys or transactional workflow boundaries are documented.
- Many scripts construct SQL using string formatting.
- `INSERT IGNORE` can hide data-quality problems.
- Several scripts can leave files and DB rows inconsistent after partial failure.

Needed work:

- Add dry-run modes before new write workflows.
- Add preflight reports for every mutating command.
- Add post-write validation queries.
- Use parameterized SQL in new code.

### Path management

- Absolute local paths are embedded in scripts and DB rows.
- Remote merlot paths are embedded in download scripts.
- Output directories are not consistently separated from source/reference data.

Needed work:

- Centralize path roots in a local config file or environment variables.
- Keep DB file paths readable but add validators for moved/missing files.
- Separate `reference_data/`, `legacy/`, `generated/`, and active source code.

### Reproducibility

- Interactive prompts make old results hard to reproduce.
- Manual DIB measurement and primary selection need review records.
- Multiple script versions exist without a clear current entry point.

Needed work:

- Add command wrappers with explicit arguments.
- Save machine-readable manifests for combine and release products.
- Store review notes separately from science DB until the schema is decided.

### Usability

- The current CLI-first workflow requires repeatedly opening generated plots.
- The spectrum viewer should become the read-only review surface before it gains
  DB write features.

Needed work:

- Add status panels for object -> observation -> reduction -> telluric ->
  combine -> DIB measurement.
- Add multi-object spectrum comparison.
- Add single-object multi-combine comparison.
- Add read-only primary-candidate inspection before DB writes.

## Immediate next tasks

1. Add `dibctl.py workflow status OBJECT_ID`.
   - Show object aliases, observations, reductions, telluric corrections,
     combines, and DIB measurement counts.
2. Add read-only preflight for observation CSV import.
3. Add read-only preflight for pipeline/reduction registration.
4. Add waveshift status check for a `pipelineID`.
5. Add telluric status check for a target `pipelineID`.
6. Add combine status check for an `objectID`.
7. Add DIB status wrapper around `DIBmeasurement_status.py`.
8. Only after those are reliable, add controlled DB write commands with
   `--dry-run` as the default behavior.
