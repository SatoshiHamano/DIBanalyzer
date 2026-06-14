# DIB Workbench DB performance notes

This note records read-only query profiling for the first spectrum viewer.
No database schema changes have been applied.

## Scope

Profile target:

- objectID: `59`
- object query: `HD147889`
- combineID: `HD_147889_o59_c1T`
- order: `42`
- repeat: `7`

Command:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 workbench_tools/profile_readonly_queries.py \
  --object-id 59 --object-query HD147889 --order 42 --repeat 7
```

## Result

| Query | Rows | Median ms | Notes |
| --- | ---: | ---: | --- |
| `objectdict_exact_lookup` | 1 | 0.160 | Uses `objectdict` primary key |
| `object_search_like_fallback` | 1 | 65.863 | Full scan and filesort |
| `object_detail` | 1 | 0.462 | Uses `object` primary key; scans aliases |
| `combine_list` | 3 | 0.270 | Scans small `combinesummary` table |
| `spectrum_orders` | 20 | 0.728 | Scans `combinedspectrum` |
| `measurements_for_order` | 1 | 4.295 | Scans `DIBmeasurement` |

## Interpretation

The current data volume is small enough that the first viewer is usable without
adding indexes.  The clear slow path is broad object search with
`LIKE '%query%'`, which scans `object` and `objectdict`.

The viewer now uses an exact alias lookup first:

```sql
select objectid
from objectdict
where registeredname = ?
order by priority
limit ?
```

This makes common searches such as `HD147889` fast without changing the
database.  If no exact alias is found, the viewer falls back to the broader
LIKE search.

## Index candidates for later

These are candidates only.  They have not been executed.

```sql
create index idx_objectdict_objectid
  on objectdict(objectid);

create index idx_combinesummary_objectID
  on combinesummary(objectID);

create index idx_combinedspectrum_combineID_order
  on combinedspectrum(combineID, echelleorder);

create index idx_DIBmeasurement_combineID_order
  on DIBmeasurement(combineID, echelleorder);
```

Expected impact:

- `objectdict(objectid)` should remove alias scans in object detail.
- `combinesummary(objectID)` should make combine listing robust if the table grows.
- `combinedspectrum(combineID, echelleorder)` should remove order lookup scans.
- `DIBmeasurement(combineID, echelleorder)` should remove measurement lookup scans.

Before applying any index, run `EXPLAIN` and timing again, then put the DDL in a
small migration SQL file instead of executing ad hoc changes.
