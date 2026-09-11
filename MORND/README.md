# Mountain Orchard Road Network Dataset (MORND)

MORND is a research-reference dataset for locating and examining mountain-orchard roads used in UAV–UGV cooperative inspection experiments. It provides stable road identifiers, ordered road points, representative geographic coordinates, and ASTER GDEM-derived elevation profiles so that researchers can locate the corresponding roads and reproduce the geometric and elevation processing used by the project.

## Scope and intended use

- Non-commercial academic research, method verification, and result reproduction.
- The road annotations were manually digitized with Google satellite imagery as a visual reference.
- The labelled map is a static research-reference image. The underlying Google GeoTIFF is not redistributed.
- Users remain responsible for complying with the applicable Google terms and attribution requirements.
- These annotations are not authoritative surveying, cadastral, safety-critical navigation, or autonomous-driving data.

## Dataset contents

- `data/mornd_roads.csv`: one row per road, including length, elevation statistics, representative coordinates, and paper road labels.
- `data/mornd_road_vertices.csv`: ordered 2 m resampled points for every road with projected coordinates, geographic coordinates, raw interpolated elevation, filled elevation, and 40 m moving-mean elevation.
- `data/mornd_road_index.csv`: compact road-ID-to-coordinate lookup table with direct Google Maps lookup URLs.
- `data/mornd_roads_processed.*`: cleaned and 2 m resampled two-dimensional Shapefile in EPSG:4542.
- `data/mornd_roads.mat`: MATLAB representation of the processed tables and geometries.
- `figures/mornd_road_id_map.jpg`: static road-identifier reference map with source attribution.
- `metadata/data_dictionary.md`: field definitions and units.
- `metadata/provenance.md`: source and processing provenance.
- `qa/qa_summary.csv` and `qa/validation_report.md`: compact data-quality evidence.
- `scripts/build_mornd.m`: reproducible build script.

## Road identifiers

- Main roads: `MORND-MAJ-001` to `MORND-MAJ-009` (`M001`–`M009` on the map).
- Field roads: `MORND-FLD-001` to `MORND-FLD-336` (`F001`–`F336` on the map).

## Coordinate and elevation reference

- Road projected coordinates: CGCS2000 / 3-degree Gauss–Krüger CM 99E, EPSG:4542, metres.
- Geographic coordinates: CGCS2000 geographic, EPSG:4490, decimal degrees.
- Elevation source: ASTER GDEM, native nominal spatial resolution 30 m.
- The working DEM file has been resampled to an approximately 10 m grid. This changes grid spacing but does not increase the effective terrain resolution beyond the 30 m source product.
- Road elevations are obtained using bilinear interpolation and then smoothed with a 40 m moving-mean window.

## Required attribution

When using the road lookup image or annotations, cite this dataset and state:

> Road annotations were manually digitized with Google satellite imagery as a visual reference. Satellite imagery © Google and imagery providers.

The Google image is not included as a reusable raster layer. Google and its imagery providers retain their respective rights in the underlying imagery.

ASTER GDEM should be cited according to the version actually used. For a clean redistributable build, the recommended source is ASTER GDEM V3 obtained directly from NASA Earthdata:

> NASA/METI/AIST/Japan Spacesystems, and U.S./Japan ASTER Science Team. ASTER Global Digital Elevation Model V003. 2018. NASA EOSDIS Land Processes DAAC. https://doi.org/10.5067/ASTER/ASTGTM.003

## Important limitation

The current build uses an ASTER GDEM working raster obtained through the Geospatial Data Cloud. Its platform terms may restrict redistribution of the raster itself. Consequently, MORND does not copy the source DEM GeoTIFF. Before depositing a fully redistributable raster, rebuild the same study-area subset from ASTER GDEM V3 downloaded directly from NASA Earthdata.(area 99.05805588°E–99.07447100°E and 24.24183848°N–24.25885937°N )
