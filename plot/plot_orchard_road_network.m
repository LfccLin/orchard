function plot_orchard_road_network(dataDir, outputDir)
%PLOT_ORCHARD_ROAD_NETWORK Create the orchard-road overview used in the paper.
% Run this file directly in MATLAB. By default, input data are read from
% <project>/dem and figure files are written to <project>/results/figures.
%
% Required toolbox: Mapping Toolbox.

close all;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
if nargin < 1 || isempty(dataDir)
    dataDir = fullfile(projectRoot, 'dem');
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(projectRoot, 'results', 'figures');
end
assert(isfolder(dataDir), 'Data directory not found: %s', dataDir);
if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Input files and reproducible selection settings
majorPath = findInput(dataDir, {'major_roads.shp', '大路.shp'});
minorPath = findInput(dataDir, {'minor_roads.shp', '小路.shp'});
demPath = findInput(dataDir, {'dem.tif', 'dem+.tif'});
imagePath = findInput(dataDir, {'orthophoto.tif', 'cgs2000.tif'});

cfg.crs = 4542;
cfg.routeLength = 500;       % m
cfg.sampleSpacing = 2;       % m
cfg.windowStride = 20;       % m
cfg.elevationSmooth = 40;    % m
cfg.eastingLimits = [5.060e5, 5.073e5];
cfg.maxAbsGradeP95 = 0.30;

%% Figure style: change legend size here
cfg.legendFontSize = 14;
cfg.axisFontSize = 15;
cfg.routeLabelFontSize = 13;
cfg.minorLineWidth = 1.15;
cfg.majorLineWidth = 2.35;
cfg.selectedLineWidth = 4.0;
cfg.minorColor = [0.57, 0.00, 0.36];
cfg.majorColor = [0.84, 0.13, 0.035];
cfg.selectedColor = [0.00, 0.31, 0.76];

%% Read geospatial data
majorShapes = shaperead(majorPath);
minorShapes = shaperead(minorPath);
[dem, Rdem] = readgeoraster(demPath, 'OutputType', 'double');
dem(dem < -10000) = NaN;
[ortho, Rortho] = readgeoraster(imagePath);
if size(ortho, 3) > 3
    ortho = ortho(:, :, 1:3);
end

roadCRS = projcrs(cfg.crs);
majorRoads = analyseShapeSet(majorShapes, Rdem, dem, roadCRS, cfg, 'Major');
minorRoads = analyseShapeSet(minorShapes, Rdem, dem, roadCRS, cfg, 'Minor');
majorCandidates = buildCandidates(majorRoads, cfg);
minorCandidates = buildCandidates(minorRoads, cfg);

majorCandidates = filterCandidates(majorCandidates, cfg);
minorCandidates = filterCandidates(minorCandidates, cfg);
assert(~isempty(majorCandidates), 'No valid 500 m main-road segment found.');
assert(numel(minorCandidates) >= 2, 'Fewer than two valid field-road segments found.');

%% Select three geometrically distinct routes
mapCentre = [mean(cfg.eastingLimits), mean(Rortho.YWorldLimits)];
mapDiagonal = hypot(diff(cfg.eastingLimits), diff(Rortho.YWorldLimits));

turnScore = robustScale([majorCandidates.totalTurn]);
reliefScore = robustScale([majorCandidates.elevRange]);
distance = arrayfun(@(r) hypot(r.centroidX-mapCentre(1), ...
    r.centroidY-mapCentre(2)), majorCandidates);
distance = distance(:)';
centrality = 1 - min(distance/(0.55*mapDiagonal), 1);
[~, indexA] = max(0.42*turnScore + 0.38*reliefScore + 0.20*centrality);
routeA = majorCandidates(indexA);

chosen = routeA;
routeB = chooseDiverse(minorCandidates, ...
    robustScale([minorCandidates.totalTurn]), chosen, mapDiagonal);
chosen(end+1) = routeB;
routeC = chooseDiverse(minorCandidates, ...
    robustScale([minorCandidates.elevRange]), chosen, mapDiagonal);

selected = [routeA, routeB, routeC];
routeIds = {'A', 'B', 'C'};
routeTypes = {'Main road', 'High curvature', 'High relief'};
for i = 1:numel(selected)
    selected(i).routeId = routeIds{i};
    selected(i).selectionType = routeTypes{i};
end

assert(all(abs([selected.length]-cfg.routeLength) < 1e-8));
assert(all(arrayfun(@(r) min(r.x) >= cfg.eastingLimits(1) && ...
    max(r.x) <= cfg.eastingLimits(2), selected)));

%% Draw the main map and three equally sized insets
figHeight = 950;
layoutBottom = 78;
layoutTop = 28;
layoutHeight = figHeight-layoutBottom-layoutTop;
panelGap = 12;
mainWidth = layoutHeight*diff(cfg.eastingLimits)/diff(Rortho.YWorldLimits);
routeCount = numel(selected);
insetHeight = (layoutHeight-(routeCount-1)*panelGap)/routeCount;
insetWidth = insetHeight;
groupWidth = mainWidth + panelGap + insetWidth;
figureMargin = 36;
figWidth = ceil(groupWidth+2*figureMargin);
fig = figure('Color', 'w', 'Position', [40, 40, figWidth, figHeight]);
mainLeft = (figWidth-groupWidth)/2;
insetLeft = mainLeft + mainWidth + panelGap;
mainPosition = [mainLeft, layoutBottom, mainWidth, layoutHeight];

mainAx = axes(fig, 'Units', 'pixels', 'Position', mainPosition, ...
    'PositionConstraint', 'innerposition');
mapshow(ortho, Rortho, 'Parent', mainAx);
hold(mainAx, 'on');
drawRoadSet(mainAx, minorShapes, cfg.minorColor, cfg.minorLineWidth);
drawRoadSet(mainAx, majorShapes, cfg.majorColor, cfg.majorLineWidth);

for i = 1:routeCount
    drawSelectedRoute(mainAx, selected(i), cfg, false);
end

axis(mainAx, 'equal');
xlim(mainAx, cfg.eastingLimits);
ylim(mainAx, Rortho.YWorldLimits);
set(mainAx, 'YDir', 'normal', 'FontName', 'Times New Roman', ...
    'FontSize', cfg.axisFontSize, 'LineWidth', 0.9, 'Layer', 'top', 'Box', 'on');
xTickValues = [cfg.eastingLimits(1):200:5.070e5, cfg.eastingLimits(2)];
xticks(mainAx, xTickValues);
xticklabels(mainAx, {'5.06', '5.062', '5.064', '5.066', ...
    '5.068', '5.07', '5.073'});
mainAx.XAxis.Exponent = 0;
mainAx.XTickLabelRotation = 0;
text(mainAx, 1.015, -0.028, '(\times10^{5})', 'Units', 'normalized', ...
    'FontName', 'Times New Roman', 'FontSize', cfg.axisFontSize, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'Clipping', 'off');
mainAx.YAxis.Exponent = 6;
xlabel(mainAx, '');
ylabel(mainAx, '');
title(mainAx, '');

hMajor = plot(mainAx, nan, nan, '-', 'Color', cfg.majorColor, ...
    'LineWidth', cfg.majorLineWidth);
hMinor = plot(mainAx, nan, nan, '-', 'Color', cfg.minorColor, ...
    'LineWidth', cfg.minorLineWidth);
hSelected = plot(mainAx, nan, nan, '-', 'Color', cfg.selectedColor, ...
    'LineWidth', cfg.selectedLineWidth);
lgd = legend(mainAx, [hMajor, hMinor, hSelected], ...
    {'Main roads', 'Field roads', 'Selected roads'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'NumColumns', 3, 'Box', 'off', 'FontName', 'Times New Roman', ...
    'FontSize', cfg.legendFontSize);
lgd.AutoUpdate = 'off';
drawnow;
lgd.Units = 'pixels';
legendPosition = lgd.Position;
legendPosition(1) = mainLeft + (mainWidth-legendPosition(3))/2;
legendPosition(2) = 14;
lgd.Position = legendPosition;
mainAx.Position = mainPosition;

drawMapAnnotations(mainAx);

routeSpans = arrayfun(@(r) max(diff([min(r.x), max(r.x)]), ...
    diff([min(r.y), max(r.y)])), selected);
insetSpan = ceil((max(routeSpans)+110)/20)*20;
for k = 1:routeCount
    insetBottom = layoutBottom + (routeCount-k)*(insetHeight+panelGap);
    ax = axes(fig, 'Units', 'pixels', ...
        'Position', [insetLeft, insetBottom, insetWidth, insetHeight], ...
        'PositionConstraint', 'innerposition');
    route = selected(k);
    [xLimits, yLimits] = squareCropLimits(route, insetSpan, Rortho);
    showOrthoCrop(ax, ortho, Rortho, xLimits, yLimits);
    hold(ax, 'on');
    drawRoadsInExtent(ax, minorShapes, xLimits, yLimits, ...
        cfg.minorColor, cfg.minorLineWidth+0.20);
    drawRoadsInExtent(ax, majorShapes, xLimits, yLimits, ...
        cfg.majorColor, cfg.majorLineWidth+0.20);
    drawSelectedRoute(ax, route, cfg, true);
    axis(ax, 'equal');
    xlim(ax, xLimits);
    ylim(ax, yLimits);
    set(ax, 'YDir', 'normal', 'XTick', [], 'YTick', [], ...
        'LineWidth', 1.0, 'Box', 'on', 'Layer', 'top');
end

pngPath = fullfile(outputDir, 'fig_orchard_road_network.png');
pdfPath = fullfile(outputDir, 'fig_orchard_road_network.pdf');
exportgraphics(fig, pngPath, 'Resolution', 300);
exportgraphics(fig, pdfPath, 'ContentType', 'image', 'Resolution', 300);
close(fig);

summary = buildSummary(selected);
writetable(summary, fullfile(outputDir, 'selected_routes.csv'));
save(fullfile(outputDir, 'selected_routes.mat'), 'selected', 'summary', 'cfg');
fprintf('Figure written to:\n  %s\n  %s\n', pngPath, pdfPath);
end

function path = findInput(dataDir, candidates)
for i = 1:numel(candidates)
    matches = dir(fullfile(dataDir, '**', candidates{i}));
    if ~isempty(matches)
        path = fullfile(matches(1).folder, matches(1).name);
        return;
    end
end
error('Missing input in %s. Expected one of: %s', ...
    dataDir, strjoin(candidates, ', '));
end

function candidates = filterCandidates(candidates, cfg)
if isempty(candidates)
    return;
end
keep = [candidates.absGradeP95] <= cfg.maxAbsGradeP95;
withinExtent = arrayfun(@(r) min(r.x) >= cfg.eastingLimits(1) && ...
    max(r.x) <= cfg.eastingLimits(2), candidates);
keep = keep(:) & withinExtent(:);
candidates = candidates(keep);
end

function roads = analyseShapeSet(shapes, Rdem, dem, roadCRS, cfg, sourceLayer)
roads = repmat(emptyRoad(), numel(shapes), 1);
for i = 1:numel(shapes)
    road = emptyRoad();
    road.sourceLayer = sourceLayer;
    road.featureId = i;
    [x, y] = longestFinitePart(shapes(i).X, shapes(i).Y);
    if numel(x) < 2
        roads(i) = road;
        continue;
    end
    duplicate = [false; hypot(diff(x), diff(y)) < 1e-6];
    x(duplicate) = [];
    y(duplicate) = [];
    s0 = [0; cumsum(hypot(diff(x), diff(y)))];
    if numel(s0) < 2 || s0(end) < 2*cfg.sampleSpacing
        roads(i) = road;
        continue;
    end
    s = (0:cfg.sampleSpacing:s0(end))';
    if s(end) < s0(end)
        s(end+1, 1) = s0(end); %#ok<AGROW>
    end
    xq = interp1(s0, x, s, 'pchip');
    yq = interp1(s0, y, s, 'pchip');
    [lat, lon] = projinv(roadCRS, xq, yq);
    [xi, yi] = geographicToIntrinsic(Rdem, lat, lon);
    zRaw = interp2(dem, xi, yi, 'linear', NaN);
    if mean(isfinite(zRaw)) < 0.90
        roads(i) = road;
        continue;
    end
    zRaw = fillmissing(zRaw, 'linear', 'EndValues', 'nearest');
    window = max(3, round(cfg.elevationSmooth/cfg.sampleSpacing));
    if mod(window, 2) == 0
        window = window+1;
    end
    z = smoothdata(zRaw, 'movmean', window);
    road = computeGeometry(road, s, xq, yq, z);
    road.zRaw = zRaw;
    road.valid = true;
    roads(i) = road;
end
end

function candidates = buildCandidates(roads, cfg)
candidates = repmat(emptyRoad(), 0, 1);
localS = (0:cfg.sampleSpacing:cfg.routeLength)';
for i = 1:numel(roads)
    road = roads(i);
    if ~road.valid || road.length < cfg.routeLength
        continue;
    end
    offsets = 0:cfg.windowStride:(road.length-cfg.routeLength);
    finalOffset = road.length-cfg.routeLength;
    if abs(offsets(end)-finalOffset) > 1e-8
        offsets(end+1) = finalOffset; %#ok<AGROW>
    end
    for j = 1:numel(offsets)
        queryS = offsets(j)+localS;
        candidate = road;
        candidate.startOffset = offsets(j);
        candidate.endOffset = offsets(j)+cfg.routeLength;
        x = interp1(road.s, road.x, queryS, 'pchip');
        y = interp1(road.s, road.y, queryS, 'pchip');
        z = interp1(road.s, road.z, queryS, 'linear');
        candidate = computeGeometry(candidate, localS, x, y, z);
        candidate.zRaw = interp1(road.s, road.zRaw, queryS, 'linear');
        candidate.valid = true;
        candidates(end+1, 1) = candidate; %#ok<AGROW>
    end
end
end

function road = chooseDiverse(pool, baseScore, chosen, mapDiagonal)
eligible = true(1, numel(pool));
used = [chosen(strcmp({chosen.sourceLayer}, 'Minor')).featureId];
if ~isempty(used)
    eligible = eligible & ~ismember([pool.featureId], used);
end
minimumDistance = inf(1, numel(pool));
for i = 1:numel(chosen)
    distance = hypot([pool.centroidX]-chosen(i).centroidX, ...
        [pool.centroidY]-chosen(i).centroidY);
    minimumDistance = min(minimumDistance, distance);
end
distanceScore = min(minimumDistance/(0.45*mapDiagonal), 1);
score = 0.74*robustScale(baseScore) + 0.26*distanceScore;
strict = eligible & minimumDistance >= 260;
if any(strict)
    eligible = strict;
end
score(~eligible) = -Inf;
[~, index] = max(score);
assert(isfinite(score(index)), 'Could not select a distinct road segment.');
road = pool(index);
end

function scaled = robustScale(values)
values = double(values(:)');
low = percentile(values, 5);
high = percentile(values, 95);
if ~isfinite(low) || ~isfinite(high) || high <= low
    scaled = 0.5*ones(size(values));
else
    scaled = min(max((values-low)/(high-low), 0), 1);
end
end

function value = percentile(values, percentage)
values = sort(values(isfinite(values)));
if isempty(values)
    value = NaN;
    return;
end
position = numel(values)*percentage/100 + 0.5;
position = min(max(position, 1), numel(values));
lower = floor(position);
upper = ceil(position);
if lower == upper
    value = values(lower);
else
    value = values(lower) + (position-lower)*(values(upper)-values(lower));
end
end

function road = computeGeometry(road, s, x, y, z)
road.s = s-s(1);
road.x = x;
road.y = y;
road.z = z;
dx = gradient(x, road.s);
dy = gradient(y, road.s);
dz = gradient(z, road.s);
road.psi = unwrap(atan2(dy, dx));
road.alpha = atan2(dz, hypot(dx, dy));
road.kappa = gradient(road.psi, road.s);
road.length = road.s(end);
road.elevRange = max(z)-min(z);
road.totalTurn = trapz(road.s, abs(road.kappa));
road.absGradeP95 = percentile(abs(tan(road.alpha)), 95);
road.centroidX = mean(x);
road.centroidY = mean(y);
end

function drawRoadSet(ax, shapes, color, lineWidth)
for i = 1:numel(shapes)
    plot(ax, shapes(i).X, shapes(i).Y, '-', ...
        'Color', color, 'LineWidth', lineWidth);
end
end

function drawSelectedRoute(ax, route, cfg, isInset)
if strcmp(route.sourceLayer, 'Major')
    routeWidth = cfg.selectedLineWidth+0.7;
else
    routeWidth = cfg.selectedLineWidth;
end
if isInset
    routeWidth = routeWidth+0.5;
end
plot(ax, route.x, route.y, '-', 'Color', 'w', 'LineWidth', routeWidth+1.6);
plot(ax, route.x, route.y, '-', 'Color', cfg.selectedColor, ...
    'LineWidth', routeWidth);
[labelX, labelY] = endpointLabelPosition(route, 32);
labelSize = cfg.routeLabelFontSize-isInset;
text(ax, labelX, labelY, route.routeId, 'Color', 'w', ...
    'FontName', 'Times New Roman', 'FontWeight', 'bold', ...
    'FontSize', labelSize, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'BackgroundColor', cfg.selectedColor, ...
    'Margin', 2.5, 'Clipping', 'on');
end

function [labelX, labelY] = endpointLabelPosition(route, offset)
direction = [route.x(2)-route.x(1), route.y(2)-route.y(1)];
direction = direction/max(norm(direction), eps);
labelX = route.x(1)-offset*direction(1);
labelY = route.y(1)-offset*direction(2);
end

function drawMapAnnotations(ax)
xl = xlim(ax);
yl = ylim(ax);
barX = xl(1)+0.055*diff(xl);
barY = yl(1)+0.055*diff(yl);
plot(ax, [barX, barX+200], [barY, barY], 'w-', ...
    'LineWidth', 5.5, 'HandleVisibility', 'off');
text(ax, barX+100, barY+0.017*diff(yl), '200 m', 'Color', 'w', ...
    'FontName', 'Times New Roman', 'FontWeight', 'bold', ...
    'FontSize', 13, 'HorizontalAlignment', 'center');
text(ax, xl(2)-0.055*diff(xl), yl(2)-0.060*diff(yl), 'N\uparrow', ...
    'Color', 'w', 'FontName', 'Times New Roman', 'FontWeight', 'bold', ...
    'FontSize', 17, 'HorizontalAlignment', 'center');
end

function [xLimits, yLimits] = squareCropLimits(route, span, R)
cx = 0.5*(min(route.x)+max(route.x));
cy = 0.5*(min(route.y)+max(route.y));
xLimits = [cx-span/2, cx+span/2];
yLimits = [cy-span/2, cy+span/2];
if xLimits(1) < R.XWorldLimits(1)
    xLimits = xLimits+(R.XWorldLimits(1)-xLimits(1));
elseif xLimits(2) > R.XWorldLimits(2)
    xLimits = xLimits-(xLimits(2)-R.XWorldLimits(2));
end
if yLimits(1) < R.YWorldLimits(1)
    yLimits = yLimits+(R.YWorldLimits(1)-yLimits(1));
elseif yLimits(2) > R.YWorldLimits(2)
    yLimits = yLimits-(yLimits(2)-R.YWorldLimits(2));
end
end

function showOrthoCrop(ax, ortho, R, xLimits, yLimits)
dx = R.CellExtentInWorldX;
dy = R.CellExtentInWorldY;
c1 = max(1, floor((xLimits(1)-R.XWorldLimits(1))/dx)+1);
c2 = min(size(ortho, 2), ceil((xLimits(2)-R.XWorldLimits(1))/dx));
r1 = max(1, floor((R.YWorldLimits(2)-yLimits(2))/dy)+1);
r2 = min(size(ortho, 1), ceil((R.YWorldLimits(2)-yLimits(1))/dy));
crop = ortho(r1:r2, c1:c2, :);
xWorld = [R.XWorldLimits(1)+(c1-1)*dx, R.XWorldLimits(1)+c2*dx];
yWorld = [R.YWorldLimits(2)-r2*dy, R.YWorldLimits(2)-(r1-1)*dy];
Rc = maprefcells(xWorld, yWorld, [size(crop, 1), size(crop, 2)], ...
    'ColumnsStartFrom', 'north');
mapshow(crop, Rc, 'Parent', ax);
end

function drawRoadsInExtent(ax, shapes, xLimits, yLimits, color, lineWidth)
for i = 1:numel(shapes)
    bounds = shapes(i).BoundingBox;
    overlaps = bounds(2,1) >= xLimits(1) && bounds(1,1) <= xLimits(2) && ...
        bounds(2,2) >= yLimits(1) && bounds(1,2) <= yLimits(2);
    if overlaps
        plot(ax, shapes(i).X, shapes(i).Y, '-', ...
            'Color', color, 'LineWidth', lineWidth);
    end
end
end

function [x, y] = longestFinitePart(X, Y)
finite = isfinite(X) & isfinite(Y);
edges = diff([false, finite, false]);
starts = find(edges == 1);
stops = find(edges == -1)-1;
if isempty(starts)
    x = [];
    y = [];
    return;
end
[~, index] = max(stops-starts+1);
x = X(starts(index):stops(index));
y = Y(starts(index):stops(index));
x = x(:);
y = y(:);
end

function summary = buildSummary(selected)
routeCount = numel(selected);
summary = table('Size', [routeCount, 9], ...
    'VariableTypes', {'string','string','string','double','double', ...
    'double','double','double','double'}, ...
    'VariableNames', {'RouteID','SelectionType','SourceLayer','FeatureID', ...
    'Length_m','ElevationRange_m','AbsGradeP95_percent', ...
    'TotalAbsoluteTurn_deg','PointCount'});
for i = 1:routeCount
    route = selected(i);
    summary.RouteID(i) = string(route.routeId);
    summary.SelectionType(i) = string(route.selectionType);
    summary.SourceLayer(i) = string(route.sourceLayer);
    summary.FeatureID(i) = route.featureId;
    summary.Length_m(i) = route.length;
    summary.ElevationRange_m(i) = route.elevRange;
    summary.AbsGradeP95_percent(i) = 100*route.absGradeP95;
    summary.TotalAbsoluteTurn_deg(i) = rad2deg(route.totalTurn);
    summary.PointCount(i) = numel(route.s);
end
end

function road = emptyRoad()
road = struct('valid', false, 'sourceLayer', '', 'featureId', 0, ...
    'routeId', '', 'selectionType', '', 's', [], 'x', [], 'y', [], ...
    'z', [], 'zRaw', [], 'psi', [], 'alpha', [], 'kappa', [], ...
    'length', 0, 'elevRange', 0, 'totalTurn', 0, ...
    'absGradeP95', Inf, 'centroidX', NaN, 'centroidY', NaN, ...
    'startOffset', 0, 'endOffset', 0);
end
