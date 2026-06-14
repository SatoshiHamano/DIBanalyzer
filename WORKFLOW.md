# DIBproject workflow map

この文書は、長く触っていなかった DIBproject の作業を再開するための
運用メモです。まずは既存スクリプトを壊さず、どの順番で何を実行していたかを
見えるようにすることを目的にしています。

DB 構造は `DB_SCHEMA_NOTES.md`、DB 登録を軸にした課題一覧は
`DB_REGISTRATION_WORKFLOW.md`、スクリプト整理方針は
`SCRIPT_ORGANIZATION.md` を参照してください。

## Big picture

大まかな流れ:

1. 天体を登録する。
2. 観測ログを DB に登録する。
3. pipeline reduction 済みデータを取得し、DB に登録する。
4. telluric 補正を行う。
5. spectrum を combine する。
6. DIB を測定する。
7. primary measurement を選ぶ。
8. summary / figure / correlation を作る。

対応する主な DB テーブル:

| Step | Main tables |
| --- | --- |
| Object registration | `object`, `objectdict` |
| Observation import | `observation` |
| Pipeline import | `datareduction`, `reducedframe` |
| Telluric correction | `telluriccorrection`, `telluricresult` |
| Combine | `combinesummary`, `combinedataset`, `combinedspectrum` |
| DIB measurement | `DIBmeasurement`, `DIBlist` |
| Literature / summary EW | `DIBEWsummary` |

## 0. Before running anything

Check DB connection:

```bash
python3 -c 'from open_mysql_project import openproject; conn,cur=openproject(); cur.execute("select database(), version()"); print(cur.fetchone()); conn.close()'
```

Check the current Git state before changing old analysis scripts:

```bash
git status --short
```

Current important caveats:

- `open_mysql_project.py` contains hard-coded DB credentials.
- There are no DB foreign key constraints; relationships are maintained by
  naming conventions and script logic.
- Many scripts read and write the DB directly.
- Several scripts are interactive and ask questions through `input()`.
- Absolute paths such as `/Users/hamano/DIB_analysis/...` and
  `/media/WD_ext/PIPELINE_DATA` are embedded in the scripts.

## 1. Register objects

Main script:

```bash
python object_add_mysql.py object_list.txt
```

Related variant:

```bash
python object_add_mysql_ver2.py object_list.txt
```

What it does:

- Reads one object name per line.
- Asks whether the list is `OBJECT` or `STANDARD`.
- Queries SIMBAD when possible.
- Inserts or updates:
  - `object`
  - `objectdict`

Useful input convention:

- `name#same` registers an alias to the previous object.

Pain points:

- Interactive prompt.
- SIMBAD/network dependency.
- DB connection is duplicated directly in this script instead of using only
  `open_mysql_project.py`.
- It is not obvious when to use `object_add_mysql.py` versus
  `object_add_mysql_ver2.py`.

## 2. Register observations

Old WODB HTML import:

```bash
python observation_mysql_update.py
```

CSV import from `observationData/pinot*`:

```bash
python observation_mysql_update_csv.py
```

What they do:

- Match target names against `objectdict.registeredname`.
- Insert or update `observation`.
- Print unregistered DIB/C2/standard targets.

Pain points:

- `observation_mysql_update.py` fetches a hard-coded WODB URL.
- `observation_mysql_update_csv.py` implicitly reads `observationData/pinot*`.
- Missing object registration blocks clean observation import.

## 3. Download / register pipeline data

Likely current script:

```bash
python datadownload_merlot_ver6.py INPUT
```

Older version:

```bash
python datadownload_merlot_ver5.py INPUT
```

What it does:

- Reads pipeline information.
- Copies reduced products from the pipeline server.
- Inserts:
  - `datareduction`
  - `reducedframe`

Important embedded paths:

- Object data root: `/Users/hamano/DIB_analysis/DIB_pipeline_dir`
- Standard data root: `/Users/hamano/DIB_analysis/Standard_pipeline_dir`
- Pipeline server root: `/media/WD_ext/PIPELINE_DATA`

Pain points:

- `ver5` and `ver6` coexist without a documented choice.
- Uses `pexpect`/scp-like interaction.
- Server, user, password, and path assumptions are embedded.
- The script both downloads files and mutates DB tables.

## 4. Telluric correction

Main script:

```bash
python telluric_auto.py TARGET_PIPELINE_ID REF_PIPELINE_ID
```

Options:

```bash
python telluric_auto.py TARGET_PIPELINE_ID REF_PIPELINE_ID --advanced
python telluric_auto.py TARGET_PIPELINE_ID REF_PIPELINE_ID --manual
python telluric_auto.py TARGET_PIPELINE_ID REF_PIPELINE_ID --fsr
```

What it does:

- Reads target and standard pipeline records from `datareduction`.
- Checks slit consistency through `reducedframe` and `observation`.
- Requires `waveshift_measure/` products under the target path.
- Writes telluric corrected spectra.
- Inserts:
  - `telluriccorrection`
  - `telluricresult`

Pain points:

- Interactive duplicate-run confirmation.
- Large logic for old and new order definitions is embedded.
- Depends on exact directory layout under the pipeline output path.
- Manual and automatic modes share the same script.

## 5. Combine spectra

Main script:

```bash
python combine_MySQL.py ID1 ID2 ...
```

The IDs can be:

- `telluricID` values from `telluriccorrection`, for telluric-corrected data.
- `pipelineID` values from `datareduction`, for uncorrected data.

What it does:

- Checks that all inputs are for the same object and mode.
- Asks interactively which frame to use for each input.
- If multiple inputs are given, asks weights interactively.
- Writes combined spectra under `/Users/hamano/DIB_analysis/DIB_pipeline_dir/...`.
- Inserts:
  - `combinesummary`
  - `combinedataset`
  - `combinedspectrum`

Outputs:

- Prints the new `combineID`.
- Saves combined FITS and a PDF inspection plot.
- Saves `telluric_spectra.pickle` next to combine products.

Pain points:

- Interactive frame and weight selection.
- Absolute output root.
- One script handles both single-dataset registration and multi-dataset combine.

## 6. DIB measurement

Main script:

```bash
python DIBanalysis.py -c COMBINE_ID -d DIBID
```

Multiple IDs:

```bash
python DIBanalysis.py -c COMBINE_ID1 COMBINE_ID2 -d DIBID1 DIBID2
```

Set primary while measuring:

```bash
python DIBanalysis.py -c COMBINE_ID -d DIBID --primary
```

Telluric transmittance threshold:

```bash
python DIBanalysis.py -c COMBINE_ID -d DIBID --telthres 0.5
```

What it does:

- Uses `DIBlist` and `combinedspectrum`.
- Writes measurement spectra/images.
- Inserts into `DIBmeasurement`.

Pain points:

- Measurement includes manual/interactive work.
- It is easy to lose track of which DIBs have been measured for which object.
- Primary selection may be done separately after inspecting results.

## 7. Check measurement status

For one object:

```bash
python DIBmeasurement_status.py OBJECT_ID
```

What it shows:

- DIBs not measured.
- DIBs measured with no primary flag.
- DIBs with multiple primary flags.
- Current primary measurement where exactly one exists.

This is one of the most useful scripts to wrap in a future `dibctl status`
command.

## 8. Set primary measurement

Single measurement:

```bash
python set_primaryflag.py MEASUREMENT_ID
```

Batch from a text file:

```bash
python set_primaryflag_list.py measurement_ids.txt
```

What it does:

- Finds the related `objectID`, `DIBID`, and `combineID`.
- Sets the selected measurement to `primaryflag=1`.
- Sets competing measurements for the same object/DIB to `primaryflag=0`.

Pain points:

- Mutates DB directly.
- No dry-run mode.
- No transaction preview.

## 9. Summary and figures

Examples:

```bash
python DIBprimary.py output.pdf -o OBJECT_ID ...
python DIBsummary_forWeakDIB.py output.pdf -o OBJECT_ID ...
python DIBcorr_NIRDIB.py output.pdf -d DIBID1 DIBID2 --norm
```

Related scripts:

- `DIBall_figure.py`
- `DIBmeasurement_figure.py`
- `DIBcorr_NIRDIB.py`
- `DIBcorr_NIRDIB_norm.py`
- `DIBcorr_forOptDIB.py`
- `correlation_dibpair.py`
- `correlation_dibpair_clustering.py`
- `correlation_ebv.py`

Pain points:

- Many figure scripts have overlapping responsibilities.
- Some currently contain local analysis adjustments.
- Output paths and object lists are often embedded or copied between scripts.

## Scripts that likely need a decision

These pairs or groups should be resolved before building a nicer interface:

| Area | Scripts | Decision needed |
| --- | --- | --- |
| Object import | `object_add_mysql.py`, `object_add_mysql_ver2.py` | Which one is current? |
| Pipeline import | `datadownload_merlot_ver5.py`, `datadownload_merlot_ver6.py` | Which one should be wrapped? |
| Observation import | `observation_mysql_update.py`, `observation_mysql_update_csv.py` | HTML WODB vs local CSV workflow |
| DIB summary | `DIBprimary.py`, `DIBsummary_*`, `DIBcorr_*` | Which outputs are still useful? |
| Telluric | `telluric_auto.py`, `telluric_manual_mysql.py`, `check_advanced_telluric.py` | Separate status/check/write operations |

## Proposed usability layer

Do not rewrite the science logic first.  Add a thin command layer that calls the
existing scripts or their functions.

Possible future commands:

```bash
python dibctl.py db status
python dibctl.py object add object_list.txt
python dibctl.py observation import-csv
python dibctl.py reduction import INPUT
python dibctl.py telluric run TARGET_PIPELINE_ID REF_PIPELINE_ID --advanced
python dibctl.py combine ID1 ID2 --frame sum --weights 1 1
python dibctl.py dib measure --combine COMBINE_ID --dib DIBID --primary
python dibctl.py dib status OBJECT_ID
python dibctl.py dib set-primary MEASUREMENT_ID
```

First targets for improvement:

1. `dibctl.py db status`
   - Verify MySQL connection.
   - Print table counts and important path roots.

2. `dibctl.py dib status OBJECT_ID`
   - Wrap `DIBmeasurement_status.py`.
   - Later add filters for unmeasured/no-primary/multiple-primary.

3. `dibctl.py combine ...`
   - Replace interactive frame/weight prompts with command-line options.

4. `dibctl.py telluric status TARGET_PIPELINE_ID`
   - Show existing telluric corrections and whether products exist on disk.

5. `dibctl.py workflow next OBJECT_ID`
   - Print the next likely action for an object based on DB state.
