function result = plot_terrain_detour_case(dataDir, routeLibraryFile, outputDir)
%PLOT_TERRAIN_DETOUR_CASE Find and plot a data-derived blocked road connection.
% The script searches endpoints in the Exp2 road library, identifies a pair
% whose direct 3-D connector violates terrain clearance, and compares two
% feasible alternatives: climb--cruise--descent and a layered-graph detour.
%
% Run directly with no inputs. Required toolbox: Mapping Toolbox.

close all;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
if nargin < 1 || isempty(dataDir)
    dataDir = fullfile(projectRoot, 'dem');
end
if nargin < 2 || isempty(routeLibraryFile)
    routeLibraryFile = fullfile(projectRoot, 'results', 'exp2', ...
        'exp2a_route_library.mat');
end
if nargin < 3 || isempty(outputDir)
    outputDir = fullfile(projectRoot, 'results', 'figures');
end
assert(isfolder(dataDir), 'Data directory not found: %s', dataDir);
assert(isfile(routeLibraryFile), ...
    'Road library not found: %s. Run run_exp2a.m first.', routeLibraryFile);
if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Reproducible settings
cfg.roadCRS = 4542;
cfg.clearance_m = 25;
cfg.directSamples = 401;
cfg.minSeparation_m = 250;
cfg.maxSeparation_m = 1100;
cfg.longitudinalStep_m = 20;
cfg.lateralStep_m = 20;
cfg.maxLateralOffset_m = 600;
cfg.maxLateralIndexChange = 1;
cfg.peakTieTolerance_m = 0.5;
cfg.topPairsToPlan = inf;
cfg.pathSampleStep_m = 5;
cfg.minimumCollisionDepth_m = 5;
cfg.minimumPeakSaving_m = 8;
cfg.roadDisplayLength_m = 180;
cfg.surfacePadding_m = 90;

%% Terrain raster and existing road tasks
demPath = findInput(dataDir, {'dem.tif', 'dem+.tif'});
[dem, Rdem] = readgeoraster(demPath, 'OutputType', 'double');
dem(dem < -10000) = NaN;
assert(ismatrix(dem) && any(isfinite(dem), 'all'), ...
    'The DEM must be a non-empty single-band elevation raster.');

loaded = load(routeLibraryFile, 'routeLibrary');
assert(isfield(loaded, 'routeLibrary') && ~isempty(loaded.routeLibrary), ...
    'routeLibrary is missing or empty.');
roads = loaded.routeLibrary(:);
roadCRS = projcrs(cfg.roadCRS);
endpoints = collectEndpoints(roads, Rdem, dem, roadCRS);
assert(min([endpoints.agl_m]) > cfg.clearance_m, ...
    'clearance_m must be below every selected road endpoint AGL.');

%% Select the endpoint pair that most clearly exposes the terrain conflict
pairs = rankBlockedPairs(endpoints, roads, Rdem, dem, roadCRS, cfg);
assert(~isempty(pairs), ...
    'No road-endpoint pair violates the requested terrain clearance.');

best = struct('score', -inf);
planCount = min(cfg.topPairsToPlan, numel(pairs));
for k = 1:planCount
    pair = pairs(k);
    pointA = endpoints(pair.endpointA).position;
    pointB = endpoints(pair.endpointB).position;
    [detour, graphInfo] = layeredDetour(pointA, pointB, Rdem, dem, ...
        roadCRS, cfg);
    if isempty(detour)
        continue;
    end

    direct = sampleStraight(pointA, pointB, cfg.directSamples);
    directGround = terrainAtXY(direct(:, 1), direct(:, 2), ...
        Rdem, dem, roadCRS);
    requiredAltitude = max(directGround+cfg.clearance_m);
    cruiseAltitude = max([requiredAltitude; pointA(3); pointB(3)]);
    raised = [pointA; pointA(1:2), cruiseAltitude; ...
        pointB(1:2), cruiseAltitude; pointB];

    raisedMetrics = routeMetrics(raised);
    detourMetrics = routeMetrics(detour);
    peakSaving = cruiseAltitude-max(detour(:, 3));
    climbSaving = raisedMetrics.climb_m-detourMetrics.climb_m;
    extraLength = detourMetrics.length3D_m-raisedMetrics.length3D_m;
    score = 0.3*pair.collisionDepth_m + 2.0*max(peakSaving, 0) + ...
        0.8*max(climbSaving, 0) - 0.02*max(extraLength, 0);

    if peakSaving >= cfg.minimumPeakSaving_m && score > best.score
        best.score = score;
        best.pair = pair;
        best.detour = detour;
        best.graphInfo = graphInfo;
        best.direct = direct;
        best.directGround = directGround;
        best.raised = raised;
        best.cruiseAltitude_m = cruiseAltitude;
        best.raisedMetrics = raisedMetrics;
        best.detourMetrics = detourMetrics;
        best.peakSaving_m = peakSaving;
        best.climbSaving_m = climbSaving;
    end
end
assert(isfield(best, 'pair'), ...
    ['Blocked endpoint pairs were found, but none produced a lower-altitude ' ...
    'detour. Increase maxLateralOffset_m or reduce minimumPeakSaving_m.']);

endpointA = endpoints(best.pair.endpointA);
endpointB = endpoints(best.pair.endpointB);
roadA = roads(endpointA.roadIndex);
roadB = roads(endpointB.roadIndex);
raisedMinimumClearance = min(best.cruiseAltitude_m-best.directGround);
detourGround = terrainAtXY(best.detour(:, 1), best.detour(:, 2), ...
    Rdem, dem, roadCRS);
detourMinimumClearance = min(best.detour(:, 3)-detourGround);
assert(best.pair.minimumDirectClearance_m < 0, ...
    'The selected direct connector does not intersect the terrain.');
assert(raisedMinimumClearance >= cfg.clearance_m-1e-6, ...
    'The raised connector violates terrain clearance.');
assert(detourMinimumClearance >= cfg.clearance_m-1e-6, ...
    'The detour violates terrain clearance.');
assert(best.peakSaving_m >= cfg.minimumPeakSaving_m, ...
    'The detour does not provide the required peak-altitude reduction.');

%% Publication figure
fig = figure('Color', 'w', 'Position', [80, 80, 1180, 760]);
ax = axes(fig);
hold(ax, 'on');

[latitudeGrid, longitudeGrid] = geographicGrid(Rdem);
[terrainX, terrainY] = projfwd(roadCRS, latitudeGrid, longitudeGrid);
allRouteX = [best.direct(:, 1); best.raised(:, 1); best.detour(:, 1)];
allRouteY = [best.direct(:, 2); best.raised(:, 2); best.detour(:, 2)];
xLimits = [min(allRouteX), max(allRouteX)] + ...
    [-cfg.surfacePadding_m, cfg.surfacePadding_m];
yLimits = [min(allRouteY), max(allRouteY)] + ...
    [-cfg.surfacePadding_m, cfg.surfacePadding_m];
surfaceMask = terrainX >= xLimits(1) & terrainX <= xLimits(2) & ...
    terrainY >= yLimits(1) & terrainY <= yLimits(2);
rowMask = any(surfaceMask, 2);
columnMask = any(surfaceMask, 1);
assert(any(rowMask) && any(columnMask), ...
    'Selected path lies outside the DEM extent.');

surf(ax, terrainX(rowMask, columnMask), terrainY(rowMask, columnMask), ...
    dem(rowMask, columnMask), dem(rowMask, columnMask), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.92, 'HandleVisibility', 'off');
colormap(ax, parula(256));
cb = colorbar(ax);
cb.Label.String = 'Terrain elevation (m)';
cb.FontName = 'Times New Roman';
cb.FontSize = 11;

hDirect = plot3(ax, best.direct(:, 1), best.direct(:, 2), ...
    best.direct(:, 3), '--', 'Color', [0.72, 0.08, 0.12], ...
    'LineWidth', 2.3, 'HandleVisibility', 'off');

% Smooth display curve for the climb--cruise--descent route.
% It keeps the same two endpoints and cruise altitude, but replaces the
% vertical climb/descent segments with gently curved transitions.
pointA = endpointA.position;
pointB = endpointB.position;
cruiseAltitude = best.cruiseAltitude_m;
transitionFraction = 0.18;   % horizontal fraction used by each curved transition
curveSamples = 60;

t = linspace(0, 1, curveSamples)';
ease = 0.5 - 0.5*cos(pi*t);  % smooth cosine easing

deltaXY = pointB(1:2) - pointA(1:2);

% Curved ascent from endpoint A.
ascXY = pointA(1:2) + (transitionFraction*t).*deltaXY;
ascZ = pointA(3) + (cruiseAltitude-pointA(3))*ease;

% Level cruise section.
cruiseStartXY = pointA(1:2) + transitionFraction*deltaXY;
cruiseEndXY = pointA(1:2) + (1-transitionFraction)*deltaXY;
cruiseSamples = 120;
tc = linspace(0, 1, cruiseSamples)';
cruiseXY = cruiseStartXY + tc.*(cruiseEndXY-cruiseStartXY);
cruiseZ = repmat(cruiseAltitude, cruiseSamples, 1);

% Curved descent into endpoint B.
descXY = pointA(1:2) + ((1-transitionFraction) + transitionFraction*t).*deltaXY;
descZ = cruiseAltitude + (pointB(3)-cruiseAltitude)*ease;

raisedDisplay = [ascXY, ascZ; ...
                 cruiseXY(2:end-1, :), cruiseZ(2:end-1); ...
                 descXY, descZ];

hRaised = plot3(ax, raisedDisplay(:, 1), raisedDisplay(:, 2), ...
    raisedDisplay(:, 3), '-', 'Color', [0.94, 0.47, 0.10], ...
    'LineWidth', 3.0, 'HandleVisibility', 'off');

hDetour = plot3(ax, best.detour(:, 1), best.detour(:, 2), ...
    best.detour(:, 3), '-', 'Color', [0.00, 0.32, 0.72], ...
    'LineWidth', 3.2, 'HandleVisibility', 'off');

plot3(ax, endpointA.position(1), endpointA.position(2), ...
    endpointA.position(3), 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0.65, 0, 0], ...
    'LineWidth', 1.2, 'HandleVisibility', 'off');
plot3(ax, endpointB.position(1), endpointB.position(2), ...
    endpointB.position(3), 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0.65, 0, 0], ...
    'LineWidth', 1.2, 'HandleVisibility', 'off');

xlabel(ax, 'Easting (m)');
ylabel(ax, 'Northing (m)');
zlabel(ax, 'Elevation (m)');
set(ax, 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'LineWidth', 0.9, 'Box', 'on', 'Layer', 'top');
axis(ax, 'tight');
view(ax, 42, 28);
grid(ax, 'on');
rotate3d(fig, 'on');

pngPath = fullfile(outputDir, 'fig_terrain_detour_case.png');
pdfPath = fullfile(outputDir, 'fig_terrain_detour_case.pdf');
exportgraphics(fig, pngPath, 'Resolution', 400);
exportgraphics(fig, pdfPath, 'ContentType', 'image', 'Resolution', 400);
% Keep the figure open for interactive rotation.

%% Auditable outputs
directMetrics = routeMetrics(best.direct);
summary = table(string(roadA.routeId), string(roadB.routeId), ...
    string(endpointA.endpointLabel), string(endpointB.endpointLabel), ...
    hypot(endpointB.position(1)-endpointA.position(1), ...
    endpointB.position(2)-endpointA.position(2)), ...
    best.pair.minimumDirectClearance_m, cfg.clearance_m, ...
    best.pair.collisionDepth_m, raisedMinimumClearance, ...
    detourMinimumClearance, directMetrics.length3D_m, ...
    best.raisedMetrics.length3D_m, ...
    best.detourMetrics.length3D_m, best.raisedMetrics.climb_m, ...
    best.detourMetrics.climb_m, best.cruiseAltitude_m, ...
    max(best.detour(:, 3)), best.peakSaving_m, ...
    'VariableNames', {'RoadA', 'RoadB', 'EndpointAType', 'EndpointBType', ...
    'HorizontalSeparation_m', 'MinimumDirectClearance_m', ...
    'RequiredClearance_m', 'DirectCollisionDepth_m', ...
    'RaisedMinimumClearance_m', 'DetourMinimumClearance_m', ...
    'DirectLength3D_m', 'RaisedLength3D_m', ...
    'DetourLength3D_m', 'RaisedClimb_m', 'DetourClimb_m', ...
    'RaisedPeakAltitude_m', 'DetourPeakAltitude_m', ...
    'PeakAltitudeReduction_m'});
writetable(summary, fullfile(outputDir, 'terrain_detour_case_summary.csv'));

result = struct();
result.config = cfg;
result.demPath = demPath;
result.routeLibraryFile = routeLibraryFile;
result.summary = summary;
result.direct = best.direct;
result.raised = best.raised;
result.detour = best.detour;
result.graphInfo = best.graphInfo;
save(fullfile(outputDir, 'terrain_detour_case.mat'), 'result');

fprintf('Selected road endpoints: %s (%s) -> %s (%s)\n', ...
    roadA.routeId, endpointA.endpointLabel, ...
    roadB.routeId, endpointB.endpointLabel);
fprintf('Direct minimum clearance: %.2f m (required %.2f m)\n', ...
    best.pair.minimumDirectClearance_m, cfg.clearance_m);
fprintf('Raised route: L3D %.2f m, climb %.2f m, peak %.2f m\n', ...
    best.raisedMetrics.length3D_m, best.raisedMetrics.climb_m, ...
    best.cruiseAltitude_m);
fprintf('Detour route: L3D %.2f m, climb %.2f m, peak %.2f m\n', ...
    best.detourMetrics.length3D_m, best.detourMetrics.climb_m, ...
    max(best.detour(:, 3)));
fprintf('Figure written to:\n  %s\n  %s\n', pngPath, pdfPath);
end

function endpoints = collectEndpoints(roads, Rdem, dem, roadCRS)
template = struct('roadIndex', 0, 'endpointIndex', 0, ...
    'endpointLabel', '', ...
    'position', zeros(1, 3), 'ground_m', NaN, 'agl_m', NaN);
endpoints = repmat(template, 2*numel(roads), 1);
count = 0;
for roadIndex = 1:numel(roads)
    points = roads(roadIndex).waypoints;
    assert(size(points, 2) == 3 && size(points, 1) >= 2, ...
        'Invalid waypoint array for route %d.', roadIndex);
    for endpointIndex = [1, size(points, 1)]
        count = count+1;
        position = points(endpointIndex, :);
        ground = terrainAtXY(position(1), position(2), Rdem, dem, roadCRS);
        endpoints(count).roadIndex = roadIndex;
        endpoints(count).endpointIndex = endpointIndex;
        endpoints(count).position = position;
        if endpointIndex == 1
            endpoints(count).endpointLabel = 'start';
        else
            endpoints(count).endpointLabel = 'end';
        end
        endpoints(count).ground_m = ground;
        endpoints(count).agl_m = position(3)-ground;
    end
end
end

function pairs = rankBlockedPairs(endpoints, roads, Rdem, dem, roadCRS, cfg)
template = struct('endpointA', 0, 'endpointB', 0, 'score', 0, ...
    'collisionDepth_m', 0, 'minimumDirectClearance_m', 0, ...
    'requiredRise_m', 0);
pairs = repmat(template, 0, 1);
for a = 1:numel(endpoints)-1
    for b = a+1:numel(endpoints)
        roadA = roads(endpoints(a).roadIndex);
        roadB = roads(endpoints(b).roadIndex);
        if endpoints(a).roadIndex == endpoints(b).roadIndex
            continue;
        end
        sameSourceRoad = strcmp(roadA.sourceLayer, roadB.sourceLayer) && ...
            roadA.featureId == roadB.featureId;
        if sameSourceRoad
            continue;
        end
        pointA = endpoints(a).position;
        pointB = endpoints(b).position;
        separation = hypot(pointB(1)-pointA(1), pointB(2)-pointA(2));
        if separation < cfg.minSeparation_m || ...
                separation > cfg.maxSeparation_m
            continue;
        end
        direct = sampleStraight(pointA, pointB, cfg.directSamples);
        ground = terrainAtXY(direct(:, 1), direct(:, 2), ...
            Rdem, dem, roadCRS);
        if any(~isfinite(ground))
            continue;
        end
        clearance = direct(:, 3)-ground;
        minimumClearance = min(clearance);
        collisionDepth = max(0, -minimumClearance);
        if collisionDepth < cfg.minimumCollisionDepth_m
            continue;
        end
        requiredRise = max(ground+cfg.clearance_m)-max(pointA(3), pointB(3));
        item = template;
        item.endpointA = a;
        item.endpointB = b;
        item.collisionDepth_m = collisionDepth;
        item.minimumDirectClearance_m = minimumClearance;
        item.requiredRise_m = requiredRise;
        item.score = collisionDepth+0.5*max(requiredRise, 0);
        pairs(end+1, 1) = item; %#ok<AGROW>
    end
end
if ~isempty(pairs)
    [~, order] = sort([pairs.score], 'descend');
    pairs = pairs(order);
end
end

function [path, info] = layeredDetour(pointA, pointB, Rdem, dem, roadCRS, cfg)
path = [];
info = struct();
delta = pointB(1:2)-pointA(1:2);
distance = hypot(delta(1), delta(2));
if distance < eps
    return;
end
forward = delta/distance;
normal = [-forward(2), forward(1)];
layerCount = max(3, ceil(distance/cfg.longitudinalStep_m));
progress = linspace(0, distance, layerCount+1)';
offsets = (-cfg.maxLateralOffset_m:cfg.lateralStep_m: ...
    cfg.maxLateralOffset_m);
[~, zeroIndex] = min(abs(offsets));
nodeCount = numel(offsets);

x = pointA(1)+progress*forward(1)+offsets*normal(1);
y = pointA(2)+progress*forward(2)+offsets*normal(2);
ground = terrainAtXY(x, y, Rdem, dem, roadCRS);
z = ground+cfg.clearance_m;
valid = isfinite(z);
valid(1, :) = false;
valid(end, :) = false;
valid(1, zeroIndex) = true;
valid(end, zeroIndex) = true;
z(1, zeroIndex) = pointA(3);
z(end, zeroIndex) = pointB(3);

peakCost = inf(layerCount+1, nodeCount);
lengthCost = inf(layerCount+1, nodeCount);
previous = zeros(layerCount+1, nodeCount, 'uint16');
peakCost(1, zeroIndex) = pointA(3);
lengthCost(1, zeroIndex) = 0;
for g = 2:layerCount+1
    for m = 1:nodeCount
        if ~valid(g, m)
            continue;
        end
        predecessorIndices = max(1, m-cfg.maxLateralIndexChange): ...
            min(nodeCount, m+cfg.maxLateralIndexChange);
        predecessorIndices = predecessorIndices(valid(g-1, predecessorIndices));
        if isempty(predecessorIndices)
            continue;
        end
        dx = x(g, m)-x(g-1, predecessorIndices);
        dy = y(g, m)-y(g-1, predecessorIndices);
        dz = z(g, m)-z(g-1, predecessorIndices);
        edgeLength = sqrt(dx.^2+dy.^2+dz.^2);
        candidatePeak = max(peakCost(g-1, predecessorIndices), z(g, m));
        candidateLength = lengthCost(g-1, predecessorIndices)+edgeLength;
        minimumPeak = min(candidatePeak);
        eligible = find(candidatePeak <= ...
            minimumPeak+cfg.peakTieTolerance_m);
        [~, eligibleIndex] = min(candidateLength(eligible));
        localIndex = eligible(eligibleIndex);
        peakCost(g, m) = candidatePeak(localIndex);
        lengthCost(g, m) = candidateLength(localIndex);
        previous(g, m) = predecessorIndices(localIndex);
    end
end
if ~isfinite(lengthCost(end, zeroIndex))
    return;
end

indices = zeros(layerCount+1, 1);
indices(end) = zeroIndex;
for g = layerCount+1:-1:2
    indices(g-1) = previous(g, indices(g));
end
linear = sub2ind(size(x), (1:layerCount+1)', indices);
coarseXY = [x(linear), y(linear)];
path = buildConstantAltitudeDetour(coarseXY, pointA, pointB, ...
    Rdem, dem, roadCRS, cfg);
if isempty(path)
    return;
end
clearance = path(:, 3)-terrainAtXY(path(:, 1), path(:, 2), ...
    Rdem, dem, roadCRS);
if any(clearance < cfg.clearance_m-1e-6)
    path = [];
    return;
end
info.layerCount = layerCount+1;
info.offsetCount = nodeCount;
info.discretePeakAltitude_m = peakCost(end, zeroIndex);
info.discreteLength3D_m = lengthCost(end, zeroIndex);
info.maximumLateralOffset_m = max(abs(offsets(indices)));
info.minimumClearance_m = min(clearance);
end

function path = buildConstantAltitudeDetour(coarseXY, pointA, pointB, ...
    Rdem, dem, roadCRS, cfg)
smoothWindow = min(7, size(coarseXY, 1));
if mod(smoothWindow, 2) == 0
    smoothWindow = smoothWindow-1;
end
if smoothWindow >= 3
    smoothedXY = smoothdata(coarseXY, 1, 'movmean', smoothWindow);
    smoothedXY([1, end], :) = coarseXY([1, end], :);
else
    smoothedXY = coarseXY;
end
segmentLengths = hypot(diff(smoothedXY(:, 1)), diff(smoothedXY(:, 2)));
progress = [0; cumsum(segmentLengths)];
if progress(end) <= 0
    path = [];
    return;
end
query = (0:cfg.pathSampleStep_m:progress(end))';
if query(end) < progress(end)
    query(end+1, 1) = progress(end);
end
x = interp1(progress, smoothedXY(:, 1), query, 'pchip');
y = interp1(progress, smoothedXY(:, 2), query, 'pchip');
ground = terrainAtXY(x, y, Rdem, dem, roadCRS);
if any(~isfinite(ground))
    path = [];
    return;
end
cruiseAltitude = max([ground+cfg.clearance_m; pointA(3); pointB(3)]);
horizontal = [x, y, repmat(cruiseAltitude, numel(x), 1)];
path = [pointA; pointA(1:2), cruiseAltitude; ...
    horizontal(2:end-1, :); pointB(1:2), cruiseAltitude; pointB];
end

function values = terrainAtXY(x, y, Rdem, dem, roadCRS)
originalSize = size(x);
[latitude, longitude] = projinv(roadCRS, x(:), y(:));
[xi, yi] = geographicToIntrinsic(Rdem, latitude, longitude);
values = interp2(dem, xi, yi, 'linear', NaN);
values = reshape(values, originalSize);
end

function points = sampleStraight(pointA, pointB, count)
fraction = linspace(0, 1, count)';
points = pointA+(pointB-pointA).*fraction;
end

function metrics = routeMetrics(points)
increments = diff(points, 1, 1);
metrics.length3D_m = sum(vecnorm(increments, 2, 2));
metrics.climb_m = sum(max(increments(:, 3), 0));
metrics.descent_m = sum(max(-increments(:, 3), 0));
end

function tail = endpointRoadTail(waypoints, endpointIndex, displayLength)
segmentLengths = [0; cumsum(vecnorm(diff(waypoints, 1, 1), 2, 2))];
if endpointIndex == 1
    keep = segmentLengths <= displayLength;
else
    keep = segmentLengths >= segmentLengths(end)-displayLength;
end
tail = waypoints(keep, :);
end

function path = findInput(rootDir, candidateNames)
for i = 1:numel(candidateNames)
    matches = dir(fullfile(rootDir, '**', candidateNames{i}));
    matches = matches(~[matches.isdir]);
    if ~isempty(matches)
        [~, order] = sort({matches.folder});
        match = matches(order(1));
        path = fullfile(match.folder, match.name);
        return;
    end
end
error('Required input not found under %s: %s', ...
    rootDir, strjoin(candidateNames, ', '));
end
