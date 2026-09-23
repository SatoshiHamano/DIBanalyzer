# Data release plan for published WINERED DIB products

この文書は、出版済み DIB 論文に対応するスペクトルや測定値をGitHub Release等で
公開できるか判断し、公開用パッケージを作るための作業メモです。

目的は、このリポジトリ全体を公開することではありません。ローカル DB、
作業途中スクリプト、未整理の中間生成物を外に出すのではなく、公開可能と
判断された reduced products だけを小さな data package として切り出します。

## Scope

想定する公開対象は、まず出版済み論文に対応するデータに限定します。

候補:

- 論文で使った reduced / wavelength-corrected / telluric-corrected spectra
- DIB 測定値テーブル
- object list
- DIB line list
- 論文 figure/table の再現に必要な補助テーブル
- データ内容を説明する README
- 論文 citation と WINERED acknowledgement

公開対象外候補:

- raw data
- 未公開プロジェクトや進行中共同研究に属するデータ
- WINERED チームのポリシー上、公開確認が必要なデータ
- MySQL dump 全体
- ローカル絶対パス、パスワード、作業ログを含むファイル
- このリポジトリ全体

## Key policy questions

公開作業の前に確認する項目:

- 出版済み論文に使った reduced spectra は外部公開してよいか。
- 公開できる処理段階はどこまでか。
  - normalized spectra
  - telluric-corrected spectra
  - combined spectra
  - measurement products
- object ごとの公開可否に例外があるか。
- FITS header に削除すべき内部情報やローカルパスが含まれるか。
- GitHub Releaseを正式な配布・引用先としてよいか。
- 将来DOIを付与する場合、論文データまたはWINERED data productのどちらとして
  記述するか。DOI付与は初回公開の前提にしない。
- ライセンスをどうするか。

### Recommended use policy after initial collaborator feedback

最終公開版は、権利者・WINEREDチームの承認を前提として、データを
`CC BY 4.0`で公開する案を第一候補とする。通常利用に共同研究・共著を必須条件と
すると、公開データとしての再利用条件が不明瞭になる。代わりに次を明示する。

- 使用したGitHub ReleaseのtagまたはURLを引用情報として示すよう求める。
- 将来dataset DOIを付与した場合は、そのDOIも推奨citationへ追加する。
- 使用した天体に対応する出版論文の引用を求める。
- 大規模再解析、系統誤差の再評価、未出版対象を使う研究では事前連絡と共同研究を
  歓迎するが、通常の再利用における法的な必須条件にはしない。
- candidate版はライセンス確定前であり、出版利用・再配布はデータ提供者への確認を
  求める。

GitHub Releaseのtagを版識別子として正式な配布・引用先にする。ライセンス、
creator、各論文との関係、推奨citationをRelease本文とREADMEの両方へ記載する。
DOI付与サービスへの登録は、永続識別子が必要になった場合の追加作業とする。

現candidate版の`MANIFEST_sanitized.csv`には`source_fits_path`としてローカル絶対
パスが含まれていた。公開manifest生成処理は、basenameとパッケージ内相対パスだけを
出力するよう修正した。GitHub Release上の既存zipは未更新なので、README改訂と
合わせてcandidate v2を作る際に差し替える。

## FITS header review

`release_inventory/candidate_release_inventory.csv` の 1000 FITS について header を
読み取り確認しました。ファイルは変更していません。

確認結果:

- files scanned: 1000
- HDU structure: primary HDU only
- unique header keys: 202
- header cards: typically 142 per FITS in the sampled files

公開前に sanitize した方がよい header:

| Key | Count | Reason |
| --- | ---: | --- |
| `FITSFILE` | 1000 | 元データ取得/検出器計算機由来の Windows path を含む |
| `OBSERVER` | 680 | `Hamano et al.` など observer 情報 |
| `WODBPI` | 680 | PI 名を含む |
| `WODBOBS` | 680 | WODB internal observation ID の可能性 |
| `WODBPROP` | 680 | WODB proposal ID の可能性 |

保持してもよさそうだが確認対象:

| Key | Count | Note |
| --- | ---: | --- |
| `WODBTHEM` | 680 | `DIB`, `DIB_add`, `DIB_profile`, `NA` |
| `PIPELINE` | 700 | pipeline version |
| `FLAT` | 1000 | calibration filename |
| `BP_MASK` | 1000 | mask filename/object-like nameを含む |
| `APSCMASK` | 1000 | aperture/scatter mask name |
| `DCLOG*` | 1000 | wavelength transform references |
| `DOPCOR*`, `VAC2AIR`, `AIR2VAC` | 1000 | wavelength/velocity correction provenance |

次の安全な作業は、元 FITS を変更せず、公開用コピーを作る段階で header から
private/internal keywords を削除する sanitizer を用意することです。

暫定 sanitizer 方針:

- remove:
  - `FITSFILE`
  - `OBSERVER`
  - `WODBPI`
  - `WODBOBS`
  - `WODBPROP`
- keep pending policy check:
  - `WODBTHEM`
  - `PIPELINE`
  - calibration/provenance keys
- add:
  - `ORIGFILE` or release manifest reference
  - `RELEASE`
  - `DOI` if a dataset DOI is assigned later
  - `REFERENC` for the paper citation

### Sanitizer script

公開用 FITS copy を作る sanitizer を用意済みです。

```text
release_tools/sanitize_fits_for_release.py
```

デフォルトでは元 FITS を変更せず、次の公開用ディレクトリにコピーを書きます。

```text
release_package/spectra/
release_package/MANIFEST_sanitized.csv
```

デフォルト削除キー:

```text
FITSFILE, OBSERVER, WODBPI, WODBOBS, WODBPROP
```

Dry-run:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 release_tools/sanitize_fits_for_release.py --dry-run --limit 3 --no-checksum
```

検証実行:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 release_tools/sanitize_fits_for_release.py --limit 45 --overwrite --no-checksum
```

検証結果:

- 45 FITS の公開用コピーを作成
- 2014年系 FITS では `FITSFILE` の削除を確認
- 2017年系 FITS では `FITSFILE`, `OBSERVER`, `WODBPI`, `WODBOBS`, `WODBPROP` の削除を確認
- `RELEASE=WINERED_DIB_PUBLISHED_DATA` を追加
- `HISTORY` に sanitizer 実行履歴と削除キーを記録

全件実行はまだしていません。次に進む場合は、reference表記を決めてから1000件
すべてを処理します。dataset DOIは未設定のまま公開可能です。

## Proposed release package

推奨する構成:

```text
winered-dib-published-data/
  README.md
  MANIFEST.csv
  CITATION.cff
  LICENSE.txt
  tables/
    objects.csv
    dib_lines.csv
    dib_measurements.csv
  spectra/
    OBJECT_NAME/
      COMBINE_ID/
        COMBINE_ID_mORDER.fits
        COMBINE_ID_mORDER.txt
  provenance/
    processing_summary.md
    source_mapping.csv
```

`MANIFEST.csv` に入れたい列:

- `object_id`
- `object_name`
- `combine_id`
- `dib_id`
- `echelle_order`
- `wavelength_min`
- `wavelength_max`
- `spectrum_file`
- `measurement_id`
- `processing_level`
- `published_reference`
- `notes`

## Repository-side work

このリポジトリ側で行う作業:

1. 論文対象 object / combine / measurement を特定する。
2. DB から公開候補 manifest を作る。
3. manifest に対応するスペクトルファイルの存在を確認する。
4. FITS header と ASCII table に公開不可情報がないか確認する。
5. 公開用ディレクトリへコピーするスクリプトを作る。
6. README と citation を作る。
7. 小さなサンプル package で内容確認する。
8. チーム確認後、版付きGitHub Releaseを作る。

## Likely source tables

DB 側の主な参照元:

- `object`
- `objectdict`
- `combinesummary`
- `combinedataset`
- `combinedspectrum`
- `DIBlist`
- `DIBmeasurement`
- `DIBEWsummary`
- `telluriccorrection`
- `telluricresult`

特に `combinedspectrum.combinefilepath` と `DIBmeasurement.DIBspecpath` は、
公開候補ファイルの実体を探す起点になります。

## Existing publication archive

`figureFactory_weakDIB/` は出版済み weak-DIB 論文用の再現性アーカイブとして
扱います。

ここにある appendix table や figure 作成スクリプトは、公開 package の
table 内容や README を確認する参考になります。ただし、このディレクトリを
そのまま公開packageに含めるのではなく、必要な成果物だけを整理して
取り出します。

## Prior shared dataset trace

過去に外部共有用にまとめた可能性があるデータセットの痕跡が見つかりました。

関連するスクリプトは `order_spec_index.xlsx` を読み、`astro_scripts_uibk_test/`
へ DIB alignment の出力を作る構成になっています。

`order_spec_index.xlsx` には次の列があります。

- `setting`
- `order`
- `star_name`
- `obs_date`
- `spec_path`
- `x_min`
- `x_max`

確認できた内容:

- 1000 rows
- 50 stars
- 50 combine datasets
- 20 echelle orders per star
- `spec_path` は `OBJECT_DIR/COMBINE_ID/COMBINE_ID_mORDER.fits` 形式

この index は、公開候補データセットの非常に強い手掛かりです。
index が指す FITS 実体は、MySQL DB の `combinedspectrum.combinefilepath`
から全件たどれることを確認済みです。

### Candidate release objects

`order_spec_index.xlsx` から復元した 50 天体:

- `BD+404220`
- `Cyg OB2 9`
- `HD 148379`
- `HD 150898`
- `HD 148184`
- `HD_147888`
- `HD 152408`
- `HD 169454`
- `Cyg OB2 10`
- `HD 144470`
- `HD 210191`
- `Hershel 36`
- `HD147889`
- `HD37742`
- `HD166937`
- `HD 214080`
- `Cl* Westerlund 1 W 33`
- `HD43384`
- `Cyg OB2 11`
- `HD 179406`
- `HD 168607`
- `HD 167264`
- `BD-16 4818`
- `Cyg OB2 8A`
- `HD 170740`
- `HD152235`
- `HD41117`
- `Cyg OB2 3`
- `HD 184915`
- `HD183143`
- `HD223385`
- `HD50064`
- `Cyg OB2 12`
- `HD20041`
- `Hen 3-1250`
- `HD 155806`
- `HD_149404`
- `Kleinmann star`
- `HD 164402`
- `HD 164353`
- `HD12953`
- `HD 185247`
- `HD 135591`
- `HD 154368`
- `HD 168625`
- `HD 151804`
- `HD 149038`
- `HD148605`
- `HD21389`
- `HD 141637`

### Next checks for the candidate dataset

1. Compare the 50-object list against `figureFactory_weakDIB/*_appendix_table.txt`.
2. Match `spec_path` combine IDs against `combinesummary.combineID`.
3. Confirm whether these spectra are exactly the already-shared dataset or only
   an analysis index derived from it.
4. Inspect headers for local paths or private metadata.

### DB match status

`order_spec_index.xlsx` 由来の 50 `combineID` は、現在の MySQL DB に全件
存在することを確認済みです。

確認結果:

- index combine IDs: 50
- matched in `combinesummary`: 50
- missing in `combinesummary`: 0
- rows in `combinedspectrum`: 1000
- distinct combine IDs in `combinedspectrum`: 50
- local spectrum files found from DB paths: 1000
- missing local spectrum files: 0

つまり、過去共有用データセットの index を起点にすれば、DB から天体情報と実 FITS
ファイルパスを全件引ける状態です。次はこの 1000 FITS を公開候補 inventory
として CSV 化し、FITS header の公開可否確認へ進むのがよいです。

### Inventory CSV

公開候補 inventory を作成済みです。

```text
release_inventory/candidate_release_inventory.csv
```

生成スクリプト:

```text
release_tools/create_release_inventory.py
```

Inventory summary:

- rows: 1000
- objects: 50
- combine IDs: 50
- echelle orders: 42-61, 20 orders per object
- release candidates: 1000
- missing files: 0
- total FITS size: 23,212,800 bytes
- paper target rows: 620 rows / 31 objects
- extra candidate rows: 380 rows / 19 objects

`paper_2206_03131_role` で、論文掲載対象と追加候補を分けています。

- `paper_2206_03131_reddened_target`: arXiv:2206.03131 の reddened target
- `extra_candidate_release_object`: 過去共有用 index に含まれる追加候補

Rigel は DIB が存在しない参照星なので、この index には含まれていません。

`publication_candidate_tags` と `publication_candidate_dois` には、他の出版済み
論文に対応する可能性が高い対象を記録します。

- `Hamano2016_CygOB2_DIB`: 2016 ApJ Cyg OB2 paper,
  DOI `10.3847/0004-637X/821/1/42`
  - object IDs: `9, 10, 11, 12, 13, 14, 15`
  - `Cyg OB2 11` は 2022 weak-DIB 論文対象ではないが、この 2016 論文対象
    として説明できる追加候補

- `Hamano2015_NIR_DIB_0p91_1p32um`: 2015 ApJ NIR DIB paper,
  DOI `10.1088/0004-637X/800/2/137`
  - source manuscript:
    `/Users/hamano/Documents/paper/DIBfirst/DIB_hamano.tex`
  - manuscript target list: 25 early-type stars plus Rigel reference
  - candidate inventory overlap: object IDs `9, 55, 75, 86, 90, 120, 123, 124, 128`
  - overlap rows: 180 rows / 9 objects
  - extra-candidate overlap: `55, 86, 120, 128`

2015 target list from the manuscript:

```text
HD2905, HD12953, HD14489, HD20041, HD21291, HD21389, HD23180,
HD24398, HD24912, HD25204, HD30614, HD36371, HD36486, HD36822,
HD37043, HD37128, HD37742, HD38771, HD41117, HD43384, HD50064,
HD190603, HD202850, HD223385, Cyg OB2 No.12
Reference: Rigel
```

主な列:

- `object_id`
- `db_object_name`
- `db_object_aliases`
- `db_sptype`
- `db_ra`
- `db_dec`
- `db_e_bv`
- `paper_2206_03131_role`
- `publication_candidate_tags`
- `publication_candidate_dois`
- `local_origin_tags`
- `local_origin_notes`
- `index_star_name`
- `combine_id`
- `echelle_order`
- `mode`
- `telluricflag`
- `combineflag`
- `obs_date`
- `db_wavelength_min`
- `db_wavelength_max`
- `local_fits_path`
- `file_exists`
- `file_size_bytes`
- `index_spec_path`
- `combine_path`

### Remaining local analysis candidates

出版済み論文の target list と直接対応しない追加候補は、`objectdict` の別名と
`object` の座標を inventory に追加して確認します。

現時点で出版論文タグなし、かつローカル解析由来として残る天体:

```text
8    BD-16 4818
25   HD 147888 / rho Oph D
26   HD 148184 / chi Oph
28   HD 149038 / mu Nor
29   HD 150898
34   HD 164353 / 67 Oph
35   HD 164402
37   HD 167264 / 15 Sgr
44   HD 184915 / kap Aql
52   HD 214080
60   HD 148605 / i Sco
65   HD 166937 / mu Sgr
133  Herschel 36
135  Kleinmann star
```

これらは `temporaly_files/NTT17c*_DIBanalysis.sh`,
`temporaly_files/pipelineID_20170730_quality.txt`,
`temporaly_files/quality_20170729.txt` などにまとまって現れるため、
`local_origin_tags=NTT17c_DIBanalysis_candidate` としました。
出版済み論文 DOI はまだ結びつけていません。

## First concrete step

最初に作るべきもの:

```bash
python dibctl.py release inventory --paper weak-dib
```

ただし、いきなり実装せず、まずは SQL で以下を確認します。

- 論文対象 object の一覧
- 各 object の primary DIB measurements
- 各 measurement が参照する `combineID`
- 各 `combineID` に対応する spectrum files
- それらのファイルがローカルに存在するか

## Open decisions

- 論文対象 object list をどこから確定するか。
  - `figureFactory_weakDIB/*_appendix_table.txt`
  - `Table1_objectlist.py`
  - DB query
  - 論文本体の table
- 公開する spectrum の処理段階をどれにするか。
- FITS と ASCII の両方を出すか。
- 測定値は DB 由来を正とするか、出版 table 由来を正とするか。
- C2 / KI / comparison 用の補助データを含めるか。

## Draft reply angle

メール返信で言えること:

- 出版済み論文に対応する reduced products を公開できるか確認したい。
- raw data や未公開共同研究データではなく、まず論文対応データに限定する。
- GitHub Releaseを版付きの正式配布元にできる。
- object list、spectra、DIB measurement table、README をまとめる方向で
  準備できる。
- 公開可否と処理段階について WINERED チーム内で確認したい。
