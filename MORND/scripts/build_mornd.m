function build_mornd(sourceRoot, outputRoot)
%BUILD_MORND Build the research-oriented Mountain Orchard Road Network Dataset.
%   BUILD_MORND(SOURCEROOT, OUTPUTROOT) reads the road shapefiles, ASTER GDEM
%   raster, and georeferenced Google satellite reference image under
%   SOURCEROOT. It writes cleaned/resampled road data, per-road summaries,
%   coordinate indexes, a labelled static reference map, and QA summaries to
%   OUTPUTROOT.
%
%   The Google image is used only to create a labelled static reference map.
%   The source GeoTIFF is not copied into the output dataset.

arguments
    sourceRoot (1,1) string
    outputRoot (1,1) string
end

sampleSpacing = 2;      % m
smoothLength = 40;      % m
roadCRS = projcrs(4542);

majorPath = findInput(sourceRoot, ["major_roads.shp", "大路.shp"], "shp");
minorPath = findInput(sourceRoot, ["minor_roads.shp", "小路.shp"], "shp");
demPath = findInput(sourceRoot, ["dem.tif", "dem+.tif"], "dem");
imagePath = findInput(sourceRoot, ["orthophoto.tif", "cgs2000.tif"], "卫星");

dataDir = fullfile(outputRoot, "data");
figureDir = fullfile(outputRoot, "figures");
qaDir = fullfile(outputRoot, "qa");
scriptDir = fullfile(outputRoot, "scripts");
for folder = [dataDir, figureDir, qaDir, scriptDir]
    if ~isfolder(folder), mkdir(folder); end
end

[dem, Rdem] = readgeoraster(demPath, "OutputType", "double");
dem(dem < -10000) = NaN;
major = shaperead(majorPath);
minor = shaperead(minorPath);

paperLabels = containers.Map( ...
    {'Major-2', 'Field-301', 'Field-14'}, {'A', 'B', 'C'});

roadRows = cell(0, 20);
vertexRows = cell(0, 13);
indexRows = cell(0, 10);
processedShapes = repmat(struct( ...
    "Geometry", "Line", "BoundingBox", [], "X", [], "Y", [], ...
    "RoadID", "", "RoadClass", "", "SrcFeatID", 0, ...
    "PaperID", "", "Length_m", 0), 0, 1);
mapRecords = repmat(struct("roadId", "", "shortId", "", "roadClass", "", ...
    "paperId", "", "x", [], "y", [], "repX", NaN, "repY", NaN), 0, 1);

sets = {major, "Major", "MAJ", "M"; minor, "Field", "FLD", "F"};
removedDuplicateCount = 0;
invalidRoadCount = 0;

for setIndex = 1:size(sets, 1)
    shapes = sets{setIndex, 1};
    roadClass = string(sets{setIndex, 2});
    idClass = string(sets{setIndex, 3});
    shortClass = string(sets{setIndex, 4});
    for featureId = 1:numel(shapes)
        [x0, y0] = longestFinitePart(shapes(featureId).X, shapes(featureId).Y);
        originalCount = numel(x0);
        if originalCount < 2
            invalidRoadCount = invalidRoadCount + 1;
            continue
        end
        duplicate = [false; hypot(diff(x0), diff(y0)) < 1e-6];
        removedDuplicateCount = removedDuplicateCount + nnz(duplicate);
        x0(duplicate) = [];
        y0(duplicate) = [];
        if numel(x0) < 2
            invalidRoadCount = invalidRoadCount + 1;
            continue
        end

        s0 = [0; cumsum(hypot(diff(x0), diff(y0)))];
        totalLength = s0(end);
        s = (0:sampleSpacing:totalLength)';
        if s(end) < totalLength - 1e-9
            s(end+1,1) = totalLength; %#ok<AGROW>
        end
        x = interp1(s0, x0, s, "pchip");
        y = interp1(s0, y0, s, "pchip");
        [lat, lon] = projinv(roadCRS, x, y);
        [xi, yi] = geographicToIntrinsic(Rdem, lat, lon);
        zRaw = interp2(dem, xi, yi, "linear", NaN);
        validFraction = mean(isfinite(zRaw));
        zFilled = fillmissing(zRaw, "linear", "EndValues", "nearest");
        window = max(3, round(smoothLength/sampleSpacing));
        if mod(window,2) == 0, window = window + 1; end
        zSmooth = smoothdata(zFilled, "movmean", window);

        roadId = string(sprintf("MORND-%s-%03d", idClass, featureId));
        shortId = string(sprintf("%s%03d", shortClass, featureId));
        paperKey = sprintf("%s-%d", roadClass, featureId);
        if isKey(paperLabels, paperKey)
            paperId = string(paperLabels(paperKey));
        else
            paperId = "";
        end

        repIndex = max(1, round(numel(s)/2));
        repLat = lat(repIndex);
        repLon = lon(repIndex);
        googleUrl = sprintf("https://www.google.com/maps/search/?api=1&query=%.8f,%.8f", ...
            repLat, repLon);
        ds3 = hypot(hypot(diff(x), diff(y)), diff(zSmooth));
        horizontal = hypot(diff(x), diff(y));
        grade = abs(diff(zSmooth)./max(horizontal, eps));
        gradeP95 = localPercentile(grade, 95) * 100;

        roadRows(end+1,:) = {roadId, shortId, roadClass, featureId, paperId, ... %#ok<AGROW>
            totalLength, originalCount, numel(s), validFraction, ...
            min(zRaw, [], "omitnan"), max(zRaw, [], "omitnan"), ...
            mean(zRaw, "omitnan"), min(zSmooth, [], "omitnan"), ...
            max(zSmooth, [], "omitnan"), max(zSmooth)-min(zSmooth), ...
            gradeP95, sum(ds3), repLon, repLat, googleUrl};

        for vertexId = 1:numel(s)
            vertexRows(end+1,:) = {roadId, roadClass, featureId, paperId, ... %#ok<AGROW>
                vertexId, s(vertexId), x(vertexId), y(vertexId), ...
                lon(vertexId), lat(vertexId), zRaw(vertexId), ...
                zFilled(vertexId), zSmooth(vertexId)};
        end

        indexRows(end+1,:) = {roadId, shortId, roadClass, featureId, paperId, ... %#ok<AGROW>
            repLon, repLat, "EPSG:4490 (CGCS2000 geographic)", googleUrl, ...
            "Manually digitized with Google satellite imagery as reference"};

        shape = struct();
        shape.Geometry = "Line";
        shape.BoundingBox = [min(x), min(y); max(x), max(y)];
        shape.X = [x(:)' NaN];
        shape.Y = [y(:)' NaN];
        shape.RoadID = char(roadId);
        shape.RoadClass = char(roadClass);
        shape.SrcFeatID = featureId;
        shape.PaperID = char(paperId);
        shape.Length_m = totalLength;
        processedShapes(end+1,1) = shape; %#ok<AGROW>

        mapRecord = struct("roadId", roadId, "shortId", shortId, ...
            "roadClass", roadClass, "paperId", paperId, "x", x, "y", y, ...
            "repX", x(repIndex), "repY", y(repIndex));
        mapRecords(end+1,1) = mapRecord; %#ok<AGROW>
    end
end

roadTable = cell2table(roadRows, "VariableNames", [ ...
    "road_id", "short_id", "road_class", "source_feature_id", "paper_road_id", ...
    "horizontal_length_m", "original_vertex_count", "resampled_vertex_count", ...
    "dem_valid_fraction", "elevation_raw_min_m", "elevation_raw_max_m", ...
    "elevation_raw_mean_m", "elevation_smooth_min_m", "elevation_smooth_max_m", ...
    "elevation_smooth_range_m", "absolute_grade_p95_percent", ...
    "three_dimensional_length_m", "representative_longitude_deg", ...
    "representative_latitude_deg", "google_maps_lookup_url"]);
roadTable.road_id = string(roadTable.road_id);
roadTable.short_id = string(roadTable.short_id);
roadTable.road_class = string(roadTable.road_class);
roadTable.paper_road_id = string(roadTable.paper_road_id);
roadTable.google_maps_lookup_url = string(roadTable.google_maps_lookup_url);
vertexTable = cell2table(vertexRows, "VariableNames", [ ...
    "road_id", "road_class", "source_feature_id", "paper_road_id", ...
    "vertex_seq", "cumulative_distance_m", "easting_cgcs2000_m", ...
    "northing_cgcs2000_m", "longitude_cgcs2000_deg", "latitude_cgcs2000_deg", ...
    "dem_elevation_raw_m", "dem_elevation_filled_m", "elevation_smoothed_m"]);
vertexTable.road_id = string(vertexTable.road_id);
vertexTable.road_class = string(vertexTable.road_class);
vertexTable.paper_road_id = string(vertexTable.paper_road_id);
indexTable = cell2table(indexRows, "VariableNames", [ ...
    "road_id", "short_id", "road_class", "source_feature_id", "paper_road_id", ...
    "representative_longitude_deg", "representative_latitude_deg", ...
    "coordinate_reference_system", "google_maps_lookup_url", "source_note"]);
indexTable.road_id = string(indexTable.road_id);
indexTable.short_id = string(indexTable.short_id);
indexTable.road_class = string(indexTable.road_class);
indexTable.paper_road_id = string(indexTable.paper_road_id);
indexTable.coordinate_reference_system = string(indexTable.coordinate_reference_system);
indexTable.google_maps_lookup_url = string(indexTable.google_maps_lookup_url);
indexTable.source_note = string(indexTable.source_note);

writetable(roadTable, fullfile(dataDir, "mornd_roads.csv"), "Encoding", "UTF-8");
writetable(vertexTable, fullfile(dataDir, "mornd_road_vertices.csv"), "Encoding", "UTF-8");
writetable(indexTable, fullfile(dataDir, "mornd_road_index.csv"), "Encoding", "UTF-8");
shapewrite(processedShapes, fullfile(dataDir, "mornd_roads_processed.shp"));
[majorFolder, majorName] = fileparts(majorPath);
majorPrj = fullfile(majorFolder, majorName + ".prj");
if isfile(majorPrj)
    copyfile(majorPrj, fullfile(dataDir, "mornd_roads_processed.prj"));
end
writelines("UTF-8", fullfile(dataDir, "mornd_roads_processed.cpg"));
save(fullfile(dataDir, "mornd_roads.mat"), "roadTable", "vertexTable", ...
    "indexTable", "processedShapes", "sampleSpacing", "smoothLength", "-v7.3");

makeLabelledMap(imagePath, mapRecords, figureDir);

qa = table(height(roadTable), sum(roadTable.road_class == "Major"), ...
    sum(roadTable.road_class == "Field"), height(vertexTable), ...
    sum(roadTable.horizontal_length_m), removedDuplicateCount, invalidRoadCount, ...
    min(roadTable.dem_valid_fraction), sum(~isfinite(vertexTable.dem_elevation_raw_m)), ...
    numel(unique(roadTable.road_id)));
qa.Properties.VariableNames = {'road_count', 'major_road_count', 'field_road_count', ...
    'vertex_count', 'total_horizontal_length_m', 'removed_consecutive_duplicates', ...
    'invalid_road_count', 'minimum_dem_valid_fraction', ...
    'missing_raw_elevation_count', 'unique_road_id_count'};
writetable(qa, fullfile(qaDir, "qa_summary.csv"), "Encoding", "UTF-8");
missingDemRoads = roadTable(roadTable.dem_valid_fraction < 1, ...
    ["road_id", "road_class", "source_feature_id", "dem_valid_fraction", ...
    "representative_longitude_deg", "representative_latitude_deg"]);
writetable(missingDemRoads, fullfile(qaDir, "roads_with_missing_raw_dem.csv"), ...
    "Encoding", "UTF-8");

sourceScript = string(mfilename("fullpath")) + ".m";
targetScript = string(fullfile(scriptDir, "build_mornd.m"));
if ~strcmpi(sourceScript, targetScript)
    copyfile(sourceScript, targetScript);
end
fprintf("MORND build complete: %d roads, %d vertices, %.3f km.\n", ...
    height(roadTable), height(vertexTable), sum(roadTable.horizontal_length_m)/1000);
end

function path = findInput(root, candidates, subfolder)
path = "";
for candidate = candidates
    probes = [fullfile(root, subfolder, candidate), fullfile(root, candidate)];
    for probe = probes
        if isfile(probe), path = probe; return; end
    end
end
error("MORND:MissingInput", "Required source file not found under %s", root);
end

function [x, y] = longestFinitePart(X, Y)
X = X(:); Y = Y(:);
finite = isfinite(X) & isfinite(Y);
changes = diff([false; finite; false]);
starts = find(changes == 1); stops = find(changes == -1)-1;
if isempty(starts), x = []; y = []; return; end
[~, idx] = max(stops-starts+1);
x = X(starts(idx):stops(idx)); y = Y(starts(idx):stops(idx));
end

function value = localPercentile(values, pct)
values = sort(values(isfinite(values)));
if isempty(values), value = NaN; return; end
position = 1 + (numel(values)-1)*pct/100;
lo = floor(position); hi = ceil(position);
if lo == hi, value = values(lo); else
    value = values(lo) + (position-lo)*(values(hi)-values(lo));
end
end

function makeLabelledMap(imagePath, roads, figureDir)
[imageData, imageRef] = readgeoraster(imagePath);
if size(imageData,3) > 3, imageData = imageData(:,:,1:3); end
fig = figure("Visible", "off", "Color", "white", "Units", "inches", ...
    "Position", [1 1 10 12]);
ax = axes(fig, "Units", "normalized", "Position", [0.08 0.08 0.84 0.86]);
mapshow(imageData, imageRef, "Parent", ax);
hold(ax, "on");
for i = 1:numel(roads)
    if roads(i).roadClass == "Major"
        color = [0.88 0.10 0.05]; width = 1.8;
    else
        color = [0.00 0.45 0.95]; width = 0.65;
    end
    plot(ax, roads(i).x, roads(i).y, "-", "Color", color, "LineWidth", width);
end
for i = 1:numel(roads)
    label = roads(i).shortId;
    if strlength(roads(i).paperId) > 0
        label = roads(i).paperId + "/" + label;
    end
    if roads(i).roadClass == "Major"
        fontSize = 7; fontWeight = "bold"; bg = [1.0 0.90 0.75];
    else
        fontSize = 4.5; fontWeight = "normal"; bg = [1.0 1.0 0.82];
    end
    text(ax, roads(i).repX, roads(i).repY, label, ...
        "FontName", "Arial", "FontSize", fontSize, "FontWeight", fontWeight, ...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle", ...
        "Color", "black", "BackgroundColor", bg, "Margin", 0.3, ...
        "Clipping", "on");
end
axis(ax, "image");
xlim(ax, imageRef.XWorldLimits); ylim(ax, imageRef.YWorldLimits);
set(ax, "YDir", "normal");
xlabel(ax, "Easting (m), CGCS2000 / EPSG:4542");
ylabel(ax, "Northing (m), CGCS2000 / EPSG:4542");
title(ax, "MORND road identifier reference map", "FontWeight", "bold");
annotation(fig, "textbox", [0.025 0.008 0.95 0.035], "String", ...
    "Satellite imagery © Google and imagery providers. Road annotations © Authors. Non-commercial research reference only.", ...
    "EdgeColor", "none", "HorizontalAlignment", "center", ...
    "FontName", "Arial", "FontSize", 11, "FontWeight", "bold");
exportgraphics(fig, fullfile(figureDir, "mornd_road_id_map.png"), ...
    "Resolution", 300, "BackgroundColor", "white");
close(fig);
end
