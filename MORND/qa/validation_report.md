# MORND validation report

The numerical QA summary is generated automatically in `qa_summary.csv`.

## Build result

- 345 roads: 9 main roads and 336 field roads.
- 36,552 ordered resampled points.
- 72.067 km total horizontal centre-line length.
- 51 consecutive duplicate or near-duplicate source vertices removed.
- 0 invalid roads and 0 duplicate road identifiers.
- 0 duplicate (`road_id`, `vertex_seq`) keys.
- 0 missing filled or smoothed elevation values.
- 67 raw DEM samples (0.183% of all points) fall immediately outside the working DEM extent; these occur on 20 boundary roads and are explicitly retained as missing in `dem_elevation_raw_m`, while `dem_elevation_filled_m` uses the nearest valid endpoint value.
- Manuscript mapping verified: A = `MORND-MAJ-002`, B = `MORND-FLD-301`, C = `MORND-FLD-014`.

## Automated checks

- unique `road_id` for every processed feature;
- expected class counts for main and field roads;
- composite vertex key (`road_id`, `vertex_seq`) is intended to be unique;
- removal of zero-length consecutive source segments;
- finite projected and geographic coordinates;
- DEM interpolation coverage by road;
- monotonic cumulative distance within each road;
- non-negative horizontal and three-dimensional lengths;
- stable mapping of manuscript Roads A, B, and C.

## Known limitations

- Road positions were manually digitized from Google satellite imagery and were not independently validated using survey control.
- ASTER GDEM has a native nominal resolution of 30 m; the approximately 10 m working grid is an interpolated grid.
- Elevation values and derived grades are not survey-grade.
- Road topology and physical connectivity have not been independently field-verified.
- The road-ID map is dense; users should open the full-resolution image and zoom to read field-road labels.
