# MORND data dictionary

## `mornd_roads.csv`

| Field | Meaning | Unit |
|---|---|---|
| `road_id` | Stable MORND road identifier | — |
| `short_id` | Compact identifier used on the reference map | — |
| `road_class` | `Major` or `Field` | — |
| `source_feature_id` | One-based record number in the original class Shapefile | — |
| `paper_road_id` | Manuscript label `A`, `B`, or `C`, otherwise blank | — |
| `horizontal_length_m` | Length in the EPSG:4542 horizontal plane | m |
| `original_vertex_count` | Vertex count before duplicate removal and resampling | count |
| `resampled_vertex_count` | Vertex count after nominal 2 m resampling | count |
| `dem_valid_fraction` | Fraction of points with a direct finite DEM interpolation | fraction |
| `elevation_raw_min_m` | Minimum directly interpolated DEM value | m |
| `elevation_raw_max_m` | Maximum directly interpolated DEM value | m |
| `elevation_raw_mean_m` | Mean directly interpolated DEM value | m |
| `elevation_smooth_min_m` | Minimum smoothed elevation | m |
| `elevation_smooth_max_m` | Maximum smoothed elevation | m |
| `elevation_smooth_range_m` | Smoothed maximum minus minimum | m |
| `absolute_grade_p95_percent` | 95th percentile of absolute smoothed longitudinal grade | % |
| `three_dimensional_length_m` | Polyline length using smoothed elevation | m |
| `representative_longitude_deg` | Longitude of the along-road midpoint | degree |
| `representative_latitude_deg` | Latitude of the along-road midpoint | degree |
| `google_maps_lookup_url` | Direct lookup URL for research convenience | URL |

## `mornd_road_vertices.csv`

The table grain is one ordered resampled point per road. The composite key is (`road_id`, `vertex_seq`).

| Field | Meaning | Unit |
|---|---|---|
| `road_id` | Stable MORND road identifier | — |
| `road_class` | `Major` or `Field` | — |
| `source_feature_id` | Original class-specific feature record | — |
| `paper_road_id` | Manuscript label, if applicable | — |
| `vertex_seq` | One-based order along the source digitizing direction | count |
| `cumulative_distance_m` | Horizontal cumulative distance from the road start | m |
| `easting_cgcs2000_m` | EPSG:4542 easting | m |
| `northing_cgcs2000_m` | EPSG:4542 northing | m |
| `longitude_cgcs2000_deg` | EPSG:4490 longitude | degree |
| `latitude_cgcs2000_deg` | EPSG:4490 latitude | degree |
| `dem_elevation_raw_m` | Bilinearly interpolated ASTER GDEM value | m |
| `dem_elevation_filled_m` | Linear/nearest-filled elevation used before smoothing | m |
| `elevation_smoothed_m` | 40 m moving-mean elevation used by the experiments | m |

The vertical datum and local elevation bias have not been independently verified by ground control. Values must be treated as ASTER GDEM-derived terrain elevations rather than survey-grade elevations.
