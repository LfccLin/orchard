# Provenance and processing record

## Road geometry

The source road centre lines were manually digitized by the authors using Google satellite imagery as a visual reference. The source image was georeferenced in CGCS2000 / 3-degree Gauss–Krüger CM 99E (EPSG:4542). The source image itself is not redistributed as a reusable raster.

Processing applied independently to each source feature:

1. select the longest finite polyline part;
2. remove consecutive points separated by less than `1e-6 m`;
3. compute cumulative horizontal arc length;
4. resample at a nominal 2 m interval using shape-preserving cubic interpolation, retaining the final endpoint;
5. convert EPSG:4542 coordinates to CGCS2000 geographic coordinates;
6. bilinearly interpolate ASTER GDEM elevation;
7. fill isolated missing values by linear interpolation with nearest end values;
8. apply a 40 m moving-mean elevation smoother;
9. calculate per-road length, relief, grade, and representative midpoint coordinates.

## Elevation

The working DEM was obtained from the Geospatial Data Cloud and identified by the authors as ASTER GDEM 30 m data. The working file has an approximately 10 m grid spacing after resampling. This is an interpolated grid spacing, not an increase in the source product's effective 30 m spatial resolution.

The source DEM raster is deliberately excluded from this package. A future fully redistributable release should rebuild the study-area raster and elevation fields from ASTER GDEM V3 downloaded directly from NASA Earthdata.

## Reference imagery

The labelled road-ID map is provided solely to help researchers visually locate roads discussed in the project. It is a static, attributed research-reference figure and is not intended to provide a reusable copy of Google imagery.
