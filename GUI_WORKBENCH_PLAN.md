# DIB Workbench plan

この文書は、DIB 測定・確認作業を軽量 GUI で回せるようにするための計画です。

現状は、CLI スクリプトで処理を実行し、PDF や画像を目で確認し、必要なら
再度 CLI を動かす運用になっています。この方式は再現性は残しやすい一方で、
天体・combine・DIB line・telluric 補正・normalization を横断して確認する
作業が重くなります。

目的は、既存の解析ロジックをすぐ置き換えることではありません。まずは
DB と FITS を読み、測定済み情報を重ねて見られる軽量 read-only viewer を作り、
手作業確認の負担を減らします。その後、review note、flag 更新、再測定ジョブの
起動へ段階的に広げます。

## Goals

- 天体、combineID、DIBID、order を素早く行き来できる。
- FITS spectrum と DB 上の DIB 測定値を同じ画面で確認できる。
- PDF を大量に開かなくても、測定状態とスペクトル形状を確認できる。
- DB の実データ、FITS パス、measurement、primary flag の関係を見える化する。
- 最初は DB を変更しない read-only GUI として安全に始める。
- 将来的には review note、測定対象キュー、再実行コマンド生成まで扱う。

## Non-goals for the first prototype

- 既存の `DIBanalysis.py` や `combine_MySQL.py` を全面的に置き換えない。
- 初期段階では GUI から DB を直接更新しない。
- 初期段階では spectrum fitting / continuum fitting を GUI 内で再実装しない。
- 初期段階では remote pipeline data download を扱わない。
- Publication-quality figure 作成は対象外。確認作業用の表示に絞る。

## Current pain points

- CLI の出力 PDF を開いて確認する作業が多い。
- 天体名、objectID、combineID、DIBID、FITS path の対応が頭に入りにくい。
- 複数スクリプトに似た目的の処理があり、どれを使うべきか迷う。
- 測定済みか、primary か、upper limit か、manual 測定かを横断的に見づらい。
- telluric 補正や combine の候補が複数ある場合に比較が重い。
- DB は外部キー制約がなく、名前規約とスクリプト側ロジックに依存している。

## Users and workflows

主な利用者は、DIBproject の DB と FITS を理解している研究者です。
まず想定する作業は次の通りです。

### Measurement review

1. object を選ぶ。
2. combineID を選ぶ。
3. DIBID または wavelength range を選ぶ。
4. spectrum を表示する。
5. DB の `DIBmeasurement` を overlay する。
6. primary measurement、EW、FWHM、center wavelength、integration range を確認する。
7. 問題がある測定を review note として記録する。

### Dataset inspection

1. 論文・候補 dataset tag を選ぶ。
2. 対象 object list を見る。
3. 各 object の combineID、order coverage、FITS existence を見る。
4. 欠損や変な path を見つける。

### Telluric / combine comparison

1. 同一 object の複数 combineID を選ぶ。
2. 同じ DIBID / wavelength range で重ねて表示する。
3. telluric flag、advanced correction、S/N を比較する。
4. primary にすべき候補を判断する。

### Multi-object comparison

1. objectID を複数指定する。
2. center wavelength と window を指定する。
3. 各 object の最初の combineID から、center wavelength を含む order を選ぶ。
4. 同じ wavelength range でスペクトルを縦積み表示する。
5. object name、objectID、combineID、order、range 内 measurement 数を確認する。

### Re-measurement preparation

1. 問題のある measurement を note する。
2. 再測定に必要な DIBID、combineID、order、range を保存する。
3. 既存 CLI で実行できるコマンド案を生成する。
4. GUI から直接 DB 更新する前に、dry-run と差分確認を挟む。

## Phased implementation

### Phase 0: DB/API inventory

目的: GUI に必要な読み取りクエリとデータ構造を固定する。

作るもの:

- `docs` または plan 上の API query list
- `workbench_tools/inspect_readonly_path.py`
- object search query
- combine list query
- DIB list query
- measurement query
- spectrum file lookup query
- index 候補メモ

確認する DB テーブル:

- `object`
- `objectdict`
- `combinesummary`
- `combinedataset`
- `combinedspectrum`
- `DIBmeasurement`
- `DIBlist`
- `telluriccorrection`
- `telluricresult`

最初に必要な query:

```sql
-- object search
select
  o.objectid,
  o.objectname,
  group_concat(od.registeredname order by od.priority separator ' | ') as aliases,
  o.sptype,
  o.ra,
  o.decli,
  o.E_BV
from object o
left join objectdict od on o.objectid = od.objectid
group by o.objectid;

-- combine list for object
select
  combineID,
  objectID,
  mode,
  telluricflag,
  combineflag,
  combinepath
from combinesummary
where objectID = ?;

-- spectrum files for combine
select
  combineID,
  echelleorder,
  combinefilepath,
  lambdamin,
  lambdamax
from combinedspectrum
where combineID = ?
order by echelleorder;

-- measurements for combine/order
select
  m.measurementID,
  m.combineID,
  m.DIBID,
  d.wavelength_air,
  d.category,
  m.echelleorder,
  m.centerlam_air,
  m.helio_velocity,
  m.EW,
  m.EWerr,
  m.FWHM,
  m.FWHMerr,
  m.integration_start,
  m.integration_end,
  m.primaryflag,
  m.autonormalizeflag,
  m.automeasurementflag,
  m.comment
from DIBmeasurement m
join DIBlist d using(DIBID)
where m.combineID = ? and m.echelleorder = ?
order by d.wavelength_air;
```

Index candidates:

```sql
create index idx_objectdict_objectid on objectdict(objectid);
create index idx_combinesummary_objectID on combinesummary(objectID);
create index idx_combinedspectrum_combineID_order on combinedspectrum(combineID, echelleorder);
create index idx_DIBmeasurement_combineID_order on DIBmeasurement(combineID, echelleorder);
create index idx_DIBmeasurement_DIBID on DIBmeasurement(DIBID);
```

これらはまだ実行しません。まず `EXPLAIN` と既存クエリ頻度を見て、
DB migration として管理できる形にします。

Phase 0 の最初の確認スクリプト:

```bash
python3 workbench_tools/inspect_readonly_path.py object-search HD147889 --limit 5
python3 workbench_tools/inspect_readonly_path.py dib-search 10780 --limit 5
python3 workbench_tools/inspect_readonly_path.py inspect 59 --order 42
PYTHONDONTWRITEBYTECODE=1 python3 workbench_tools/profile_readonly_queries.py --object-id 59 --object-query HD147889 --order 42 --repeat 7
```

このスクリプトは `SELECT` のみを実行します。GUI の最小 API に必要な
`object -> combine -> order -> FITS path -> DIBmeasurement` の経路確認を
目的としており、DB や FITS には書き込みません。

DB 性能確認の結果は `WORKBENCH_DB_PERFORMANCE.md` にまとめます。
現時点では index は追加せず、`objectdict.registeredname` の完全一致検索を
先に試し、見つからない場合だけ broad LIKE search に戻る方針にします。

DB 読み取り処理は `workbench_tools/db_readonly.py` に集約します。
CLI inspector、query profiler、spectrum viewer はこのモジュールを共通利用し、
同じ SQL と同じ返却形式を見るようにします。

### Phase 1: Read-only spectrum viewer

目的: PDF を開かず、DB 上の測定値と FITS spectrum を同じ画面で確認する。

最初の実装:

- `workbench_tools/spectrum_viewer.py`
- `workbench_tools/spectrum_io.py`
- Python 標準 HTTP server
- `astropy.io.fits` による FITS 読み込み
- Canvas による spectrum 描画
- DB は `SELECT` のみ
- FITS は読み込みのみ

起動例:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 workbench_tools/spectrum_viewer.py --host 127.0.0.1 --port 8766
```

最小機能:

- object search
- object detail
- combine list
- order list
- spectrum plot
- DIBmeasurement overlay
- DIB line marker
- measurement detail panel
- measurement click zoom
- preserve wavelength range across combine/order changes where possible

画面:

- 左 sidebar: object search / dataset filter
- 中央: spectrum plot
- 右 panel: combine metadata / measurement table
- 下部: selected measurement detail

Plot 表示:

- wavelength vs flux
- DIB rest wavelength vertical marker
- measured center marker
- integration range shaded region
- primary measurement highlight
- telluric-heavy region marker if available

操作性メモ:

- object を切り替えたときは wavelength range を auto に戻す。
- combine を切り替えたときは、現在の wavelength range を維持する。
- order を切り替えたときも、手入力または measurement zoom の range を維持する。
- combine 変更時に現在 order が存在しない場合は、表示中 wavelength range の中心を
  coverage に含む order を優先する。
- measurement をクリックしたときは、integration range に余白を付けて zoom する。
- measurement list には状態バッジを出す。
- 選択中 measurement の詳細は Selection panel に出す。
- `EW=0`、`FWHM=0`、low S/N、comment ありなどは flagged として数える。

安全条件:

- DB write はしない。
- FITS write はしない。
- 表示用 downsample は memory 上だけで行う。
- 読み込んだ FITS path と DB query を画面に表示できるようにする。

技術候補:

- Backend: FastAPI
- Frontend: React + Plotly
- FITS reading: `astropy.io.fits`
- Data transport: JSON, large spectrum は downsample 済み配列

FITS 読み込み、線形 wavelength 復元、finite filtering、波長範囲切り出し、
downsample、表示用 robust flux limit は `workbench_tools/spectrum_io.py` に
集約します。viewer、CLI preview、将来の batch sanity check はこの共通処理を
使います。

Python だけで早く試す代替案:

- Streamlit
- Dash

長期的には、操作性と状態管理のため Web frontend を推奨します。

### Phase 2: Review notes

目的: 目視確認の結果を DB 本体に直接書かずに蓄積する。

作るもの:

- `review_notes` 用の lightweight storage
- `workbench_tools/review_notes.py`
- `workbench_review_notes.json`
- viewer API:
  - `GET /api/review-notes?measurement_ids=...`
  - `POST /api/review-note`

初期 status:

- `ok`
- `check`
- `bad_continuum`
- `telluric`
- `remeasure`
- `ignore`

保存キーは `measurementID` とする。保存内容には status、note、object/combine/order/DIB
の context、created/updated timestamp を入れる。DB には書き込まない。
- 初期は SQLite または CSV/JSONL
- 後で MySQL table に移行可能な schema

記録する情報:

- objectID
- combineID
- echelleorder
- DIBID
- measurementID
- note type
- note text
- reviewer
- created_at
- status

note type 例:

- `bad_continuum`
- `telluric_residual`
- `stellar_blend`
- `wrong_primary`
- `needs_remeasure`
- `upper_limit_check`
- `good`

この段階でも、既存の `DIBmeasurement` は変更しません。

### Phase 3: Comparison views

目的: 同じ DIB を複数天体・複数 combine で比較する。

最初の multi-object comparison 実装:

- viewer API:
  - `GET /api/compare?object_ids=59,66,9&wavelength=13175.9&window=25`
  - `GET /api/compare-combines?object_id=59&wavelength=13175.9&window=25`
- 左 sidebar の Compare panel
- objectID リスト、center wavelength、window 入力
- single objectID 入力
- 各 object の最初の combineID を自動選択
- single object mode では object の全 combineID を比較
- center wavelength を含む echelle order を自動選択
- 同一 wavelength range の縦積み Canvas plot
- DIB center marker と integration range overlay

初期制約:

- combineID は自動選択のみ。
- objectID 指定のみ。別名検索からの複数選択は後で追加する。
- DB と FITS は read-only。

機能:

- same DIBID across objects
- same object across combineIDs
- primary vs non-primary comparison
- advanced telluric vs normal telluric comparison
- normalized / raw-ish display mode where available

表示:

- velocity-space plot
- stacked spectra
- EW vs E(B-V)
- S/N and flags table

この段階で、既存の `DIBcorr_*`, `DIBsummary_*`, `correlation_*` の一部結果を
GUI から参照できるようにする。

### Phase 4: Controlled DB updates

目的: GUI から review 結果にもとづく最小限の DB 更新を安全に行う。

対象候補:

- `primaryflag`
- measurement `comment`
- review status table
- remeasurement queue table

必須条件:

- すべて dry-run first
- 更新前後の diff 表示
- 更新ログを保存
- 可能なら transaction 使用
- undo 用 SQL を保存

この段階に入る前に DB schema と index を整理します。

### Phase 5: Job launcher

目的: GUI で選んだ対象に対し、既存 CLI の実行コマンドを安全に作る。

最初は実行せず、コマンド生成だけ:

```bash
python DIBanalysis.py -d DIBID -c COMBINEID
python DIBsummary_forWeakDIB.py ...
python combine_MySQL.py ...
```

次に、実行キュー:

- pending
- running
- succeeded
- failed
- skipped

GUI 側から直接古いスクリプトを呼ぶ場合は、stdout/stderr、作成ファイル、
DB 変更を必ず記録します。

## Proposed architecture

```text
dib_workbench/
  backend/
    app.py
    db.py
    models.py
    spectra.py
    queries.py
  frontend/
    package.json
    src/
      App.tsx
      components/
      api/
  README.md
```

Backend responsibilities:

- DB connection
- query endpoints
- FITS loading
- spectrum downsampling
- path validation
- read-only API

Frontend responsibilities:

- object search
- navigation state
- plot interactions
- measurement overlays
- note UI

API sketch:

```text
GET /api/health
GET /api/objects?query=HD147889
GET /api/objects/{object_id}
GET /api/objects/{object_id}/combines
GET /api/combines/{combine_id}/orders
GET /api/combines/{combine_id}/orders/{order}/spectrum
GET /api/combines/{combine_id}/orders/{order}/measurements
GET /api/dibs?query=10780
GET /api/dibs/{dib_id}
```

Possible later write APIs:

```text
POST /api/review-notes
POST /api/remeasurement-queue
POST /api/measurements/{measurement_id}/primaryflag/dry-run
POST /api/measurements/{measurement_id}/primaryflag/apply
```

## Spectrum loading strategy

FITS files are small enough for local reads, but GUI interaction should still be
responsive.

Rules:

- Read only selected order file.
- Cache recently opened spectra in memory.
- Downsample for browser display when points exceed a threshold.
- Preserve full-resolution data for local zoom when needed.
- Return wavelength, flux, and header summary.

Open questions:

- FITS data layout is not yet standardized in this plan.
- Need to confirm whether wavelength is encoded in header WCS, table columns,
  or helper functions in `Spec1Dtools.py`.
- Need to support both WIDE and HIRES modes later.

## DB optimization plan

GUI will surface slow queries quickly. Before adding write features, DB indexes
should be reviewed.

Initial checks:

```sql
show index from DIBmeasurement;
show index from combinedspectrum;
show index from combinesummary;
show index from objectdict;
explain <gui query>;
```

Likely index candidates:

- `objectdict(objectid)`
- `objectdict(registeredname)`
- `combinesummary(objectID)`
- `combinedspectrum(combineID, echelleorder)`
- `DIBmeasurement(combineID, echelleorder)`
- `DIBmeasurement(DIBID)`
- `DIBmeasurement(primaryflag)`

Index changes should be stored as scripts, not typed manually into the DB.

Proposed directory:

```text
db_migrations/
  001_add_workbench_read_indexes.sql
```

Policy:

- First write `EXPLAIN` before/after notes.
- Keep migration idempotent where MySQL version allows.
- Do not add foreign key constraints until legacy data consistency is audited.

## First prototype milestone

The first useful milestone is:

```text
objectID -> combineID -> order -> spectrum + DIBmeasurement overlay
```

Acceptance criteria:

- Start local app with one command.
- Search object by objectID, primary name, or alias.
- Select a combineID.
- Select an echelle order.
- Show FITS spectrum in browser.
- Overlay all measurements for that combine/order.
- Highlight primary measurements.
- Show measurement detail when clicking a marker or table row.
- No DB writes.
- No FITS writes.

Suggested test objects:

- `59` / HD 147889
- `66` / HD 183143
- `9` / Cyg OB2 No.12
- `135` / Kleinmann
- `138` / Rigel

## Risks

- FITS wavelength loading may differ by product generation path.
- Some DB paths may point to files on external disks or missing mounts.
- Existing scripts may encode important correction logic outside DB fields.
- GUI can accidentally become a second analysis implementation if scope is not
  controlled.
- Write features can corrupt analysis state unless dry-run and logging are
  mandatory.

## Decisions to make

- FastAPI + React or Streamlit first.
- Where to store review notes.
- Whether GUI project should live inside this repository or a sibling repo.
- Which indexes are safe to add first.
- Which DB writes, if any, should be allowed from GUI.
- Whether remeasurement should be controlled by queue files or direct subprocess.

## Recommended next steps

1. Confirm FITS loading for 3 representative files using `Spec1Dtools.py` or
   direct `astropy.io.fits`.
2. Write read-only backend endpoints for object, combine, order, measurement.
3. Build a minimal spectrum plot page.
4. Add measurement overlay and primary highlight.
5. Add object alias and coordinate display.
6. Add review note storage as a separate local file or SQLite DB.
7. Only after this is useful, discuss DB write operations.
