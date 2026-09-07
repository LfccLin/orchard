function run_exp1()
%RUN_EXP1 Experiment 1: footprint-constrained waypoint generation.
% Run this file directly in MATLAB. It compares three methods on the same
% three 500 m orchard roads and writes publication figures plus numerical
% results to <project>/results/exp1.
%
% Method (a): terrain-normal camera pose, fixed image orientation, and the
%             largest globally feasible uniform candidate interval.
% Method (b): terrain-normal camera pose, road-aligned image orientation,
%             and the largest globally feasible uniform candidate interval.
% Method (c): the same road-aligned camera pose as (b), with the proposed
%             farthest-feasible adaptive observation-point selection.
%
% No empirical turn threshold is used. All methods share the camera, the
% horizontally transverse road surface, the 2 m candidate set, the overlap
% threshold, and the full-coverage rule.

close all;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
routeFile = fullfile(projectRoot, 'results', 'figures', ...
    'selected_routes.mat');
outputDir = fullfile(projectRoot, 'results', 'exp1');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

if ~isfile(routeFile)
    plotDir = fullfile(projectRoot, 'plot');
    addpath(plotDir);
    plot_orchard_road_network(fullfile(projectRoot, 'dem'), ...
        fullfile(projectRoot, 'results', 'figures'));
end

inputData = load(routeFile, 'selected', 'cfg');
selected = inputData.selected(1:3);
selectionCfg = inputData.cfg;
assert(numel(selected) == 3, 'Expected roads A, B, and C.');

%% Parameters fixed in the manuscript
param.routeLength = 500;              % m
param.centerlineSpacing = 2;          % m
param.mainRoadWidth = 4.5;            % m
param.fieldRoadWidth = 3.5;           % m
param.cameraDistance = 30;            % m
param.fovAlong = deg2rad(73.68);      % rad
param.fovAcross = deg2rad(53.08);     % rad
param.minimumOverlap = 0.80;
param.minimumRoadCoverage = 1.00;
param.footprintAlong = 2*param.cameraDistance* ...
    tan(param.fovAlong/2);
param.footprintAcross = 2*param.cameraDistance* ...
    tan(param.fovAcross/2);
param.maximumSpacing = (1-param.minimumOverlap)*param.footprintAlong;

assert(abs(selectionCfg.routeLength-param.routeLength) < 1e-10);
assert(abs(selectionCfg.sampleSpacing-param.centerlineSpacing) < 1e-10);
assert(abs(param.footprintAlong-44.95) < 0.03);
assert(abs(param.maximumSpacing-8.99) < 0.01);
assert(param.footprintAcross >= param.mainRoadWidth);

%% One road surface shared by all methods
% At every centreline station, the road-width direction is parallel to the
% global XY plane. Elevation therefore varies only along the centreline.
roadModelCells = cell(1, numel(selected));
for roadIndex = 1:numel(selected)
    if strcmp(char(selected(roadIndex).routeId), 'A')
        roadWidth = param.mainRoadWidth;
    else
        roadWidth = param.fieldRoadWidth;
    end
    roadModelCells{roadIndex} = buildRoadModel( ...
        selected(roadIndex), roadWidth, param);
end
roadModels = [roadModelCells{:}];

%% Three comparison methods
methods(1) = struct('id', 'a', 'name', 'Fixed orientation', ...
    'orientationMode', 'fixed', 'samplingMode', 'uniform', ...
    'color', [0.000, 0.447, 0.741]);
methods(2) = struct('id', 'b', 'name', 'Road aligned', ...
    'orientationMode', 'road', 'samplingMode', 'uniform', ...
    'color', [0.850, 0.325, 0.098]);
methods(3) = struct('id', 'c', 'name', 'Proposed adaptive', ...
    'orientationMode', 'road', 'samplingMode', 'adaptive', ...
    'color', [0.180, 0.620, 0.360]);

resultCells = cell(numel(methods), numel(roadModels));
summaryRows = cell(0, 12);
for roadIndex = 1:numel(roadModels)
    model = roadModels(roadIndex);

    fixedPoseBank = buildPoseBank(model, methods(1), param);
    alignedPoseBank = buildPoseBank(model, methods(2), param);
    fixedCoverageBank = buildCoverageBank(model, fixedPoseBank, param);
    alignedCoverageBank = buildCoverageBank(model, alignedPoseBank, param);
    [localMaximumFeasibleSpacing, localFarthestFeasibleIndex] = ...
        computeLocalMaximumFeasibleSpacing(model, alignedCoverageBank, param);

    for methodIndex = 1:numel(methods)
        method = methods(methodIndex);
        if strcmp(method.orientationMode, 'fixed')
            poseBank = fixedPoseBank;
            coverageBank = fixedCoverageBank;
        else
            poseBank = alignedPoseBank;
            coverageBank = alignedCoverageBank;
        end

        if strcmp(method.samplingMode, 'uniform')
            [selectedIndices, selectionFeasible] = selectUniformCandidates( ...
                model, coverageBank, param);
        elseif strcmp(method.samplingMode, 'adaptive')
            selectedIndices = selectAdaptiveCandidates( ...
                model, coverageBank, param);
            selectionFeasible = true;
        else
            error('Unknown sampling mode: %s', method.samplingMode);
        end

        result = assembleResult(model, poseBank, coverageBank, ...
            selectedIndices, param);
        result.constraintFeasible = selectionFeasible;
        % These candidate-wise quantities use the road-aligned footprint
        % bank shared by methods (b) and (c). They are saved for Exp. 1c
        % only and do not change waypoint selection in this experiment.
        result.localMaximumFeasibleSpacing = localMaximumFeasibleSpacing;
        result.localFarthestFeasibleIndex = localFarthestFeasibleIndex;
        validateResult(result, model, method, param);
        resultCells{methodIndex, roadIndex} = result;

        summaryRows(end+1, :) = {string(method.id), ...
            string(model.routeId), numel(result.selectedIndices), ...
            result.constraintFeasible, mean(result.intervalDistance), ...
            min(result.intervalDistance), ...
            max(result.intervalDistance), 100*min(result.adjacentOverlap), ...
            100*min(result.intervalCoverage), ...
            100*result.totalRoadCoverage, ...
            100*mean(result.footprintUtilization), ...
            mean(vecnorm(result.waypoints-result.groundPoints, 2, 2))}; ...
            %#ok<AGROW>
    end
end

results = repmat(struct('methodId', '', 'methodName', '', 'roads', []), ...
    1, numel(methods));
for methodIndex = 1:numel(methods)
    results(methodIndex).methodId = methods(methodIndex).id;
    results(methodIndex).methodName = methods(methodIndex).name;
    results(methodIndex).roads = [resultCells{methodIndex, :}];
end

summary = cell2table(summaryRows, 'VariableNames', ...
    {'Method', 'Road', 'WaypointCount', 'ConstraintFeasible', 'MeanSpacing_m', ...
    'MinSpacing_m', 'MaxSpacing_m', 'MinAdjacentOverlap_pct', ...
    'MinIntervalCoverage_pct', 'TotalRoadCoverage_pct', ...
    'MeanFootprintUtilization_pct', 'MeanCameraGroundDistance_m'});

comparisonCells = cell(1, numel(roadModels));
for roadIndex = 1:numel(roadModels)
    roadResults = [resultCells{:, roadIndex}];
    comparisonCells{roadIndex} = findComparisonLocation( ...
        roadModels(roadIndex), roadResults, methods, param);
end
comparisons = [comparisonCells{:}];

writetable(summary, fullfile(outputDir, 'exp1_summary.csv'));
save(fullfile(outputDir, 'exp1_results.mat'), ...
    'roadModels', 'methods', 'param', 'results', 'summary', 'comparisons');

%% Publication figures: one road per figure, three methods in one row
for roadIndex = 1:numel(roadModels)
    roadResults = [resultCells{:, roadIndex}];
    drawRoadComparisonFigure(roadModels(roadIndex), roadResults, ...
        methods, comparisons(roadIndex), outputDir);
end

fprintf('Experiment 1 completed. Results written to:\n  %s\n', outputDir);
disp(summary(:, {'Method', 'Road', 'WaypointCount', ...
    'MinAdjacentOverlap_pct', 'MinIntervalCoverage_pct', ...
    'TotalRoadCoverage_pct'}));
end

function model = buildRoadModel(route, roadWidth, param)
s = route.s(:);
x = route.x(:);
y = route.y(:);
z = route.z(:);
expectedCount = round(param.routeLength/param.centerlineSpacing)+1;
assert(numel(s) == expectedCount);
assert(all(isfinite([s; x; y; z])));
assert(abs(s(end)-s(1)-param.routeLength) < 1e-8);

% A rigid horizontal transform makes each road compact in the figure while
% preserving every distance, angle, slope, and coverage calculation.
xy = [x, y];
origin = mean(xy, 1);
covariance = cov(xy);
[eigenvectors, eigenvalues] = eig(covariance);
[~, principalIndex] = max(diag(eigenvalues));
axisAlong = eigenvectors(:, principalIndex);
if dot(axisAlong, (xy(end, :)-xy(1, :))') < 0
    axisAlong = -axisAlong;
end
axisAcross = [-axisAlong(2); axisAlong(1)];
localXY = [(xy-origin)*axisAlong, (xy-origin)*axisAcross];
localZ = z-min(z);
centreline = [localXY, localZ];

dx = gradient(x, s);
dy = gradient(y, s);
dz = gradient(z, s);
heading = unwrap(atan2(dy, dx));
horizontalRate = hypot(dx, dy);
beta = atan2(dz, horizontalRate);      % longitudinal slope, manuscript beta
alpha = zeros(size(beta));            % transverse slope, manuscript alpha

% Manuscript local frame:
% t = [cos(beta)cos(psi), cos(beta)sin(psi), sin(beta)]
% r = [-sin(psi), cos(psi), 0], n = t x r.
tangentWorld = [cos(beta).*cos(heading), ...
    cos(beta).*sin(heading), sin(beta)];
lateralWorld = [-sin(heading), cos(heading), zeros(size(heading))];
tangent = [tangentWorld(:, 1:2)*axisAlong, ...
    tangentWorld(:, 1:2)*axisAcross, tangentWorld(:, 3)];
lateral = [lateralWorld(:, 1:2)*axisAlong, ...
    lateralWorld(:, 1:2)*axisAcross, lateralWorld(:, 3)];
tangent = tangent./vecnorm(tangent, 2, 2);
lateral = lateral./vecnorm(lateral, 2, 2);
normal = cross(tangent, lateral, 2);
normal = normal./vecnorm(normal, 2, 2);
qAlong = tangent;
qAcross = lateral;

% The seven lateral vertices are the existing road-surface discretization.
% The normalized transverse vector keeps the physical road width exact.
lateralOffsets = linspace(-roadWidth/2, roadWidth/2, 7);
surfaceX = centreline(:, 1)+lateral(:, 1)*lateralOffsets;
surfaceY = centreline(:, 2)+lateral(:, 2)*lateralOffsets;
surfaceZ = centreline(:, 3)+lateral(:, 3)*lateralOffsets;

segmentDistance = vecnorm(diff(centreline, 1, 1), 2, 2);
[triangleV1, triangleV2, triangleV3, triangleAreas, ...
    triangleLongitudinalIndex] = roadSurfaceTriangles( ...
    surfaceX, surfaceY, surfaceZ);

model.routeId = char(route.routeId);
model.selectionType = char(route.selectionType);
model.width = roadWidth;
model.s = s;
model.worldCoordinates = [x, y, z];
model.localCoordinates = centreline;
model.heading = heading;
model.longitudinalSlope = beta;
model.crossSlope = alpha;
model.qAlong = qAlong;
model.qAcross = qAcross;
model.tangent = tangent;
model.lateral = lateral;
model.normal = normal;
model.segmentDistance = segmentDistance;
model.cumulativeDistance = [0; cumsum(segmentDistance)];
model.surfaceX = surfaceX;
model.surfaceY = surfaceY;
model.surfaceZ = surfaceZ;
model.surfaceElevation = surfaceZ+min(z);
model.triangleV1 = triangleV1;
model.triangleV2 = triangleV2;
model.triangleV3 = triangleV3;
model.triangleAreas = triangleAreas;
model.triangleLongitudinalIndex = triangleLongitudinalIndex;
model.horizontalOrigin = origin;
model.horizontalAxes = [axisAlong, axisAcross];

assert(max(abs(alpha)) < 1e-12);
assert(max(abs(lateral(:, 3))) < 1e-12);
assert(max(abs(sum(tangent.*lateral, 2))) < 1e-12);
assert(max(abs(vecnorm(normal, 2, 2)-1)) < 1e-12);
assert(all(normal(:, 3) > 0));
assert(abs(sum(triangleAreas)-sumRoadSurfaceArea(model)) < 1e-9);
end

function [v1, v2, v3, areas, longitudinalIndex] = ...
    roadSurfaceTriangles(surfaceX, surfaceY, surfaceZ)
longitudinalCount = size(surfaceX, 1)-1;
lateralCount = size(surfaceX, 2)-1;
triangleCount = 2*longitudinalCount*lateralCount;
v1 = zeros(triangleCount, 3);
v2 = zeros(triangleCount, 3);
v3 = zeros(triangleCount, 3);
areas = zeros(triangleCount, 1);
longitudinalIndex = zeros(triangleCount, 1);
triangleIndex = 0;
for longitudinal = 1:longitudinalCount
    for lateral = 1:lateralCount
        vertices = [surfaceX(longitudinal, lateral), ...
            surfaceY(longitudinal, lateral), surfaceZ(longitudinal, lateral); ...
            surfaceX(longitudinal+1, lateral), ...
            surfaceY(longitudinal+1, lateral), surfaceZ(longitudinal+1, lateral); ...
            surfaceX(longitudinal+1, lateral+1), ...
            surfaceY(longitudinal+1, lateral+1), ...
            surfaceZ(longitudinal+1, lateral+1); ...
            surfaceX(longitudinal, lateral+1), ...
            surfaceY(longitudinal, lateral+1), surfaceZ(longitudinal, lateral+1)];

        triangleIndex = triangleIndex+1;
        v1(triangleIndex, :) = vertices(1, :);
        v2(triangleIndex, :) = vertices(2, :);
        v3(triangleIndex, :) = vertices(3, :);
        areas(triangleIndex) = triangleArea(v1(triangleIndex, :), ...
            v2(triangleIndex, :), v3(triangleIndex, :));
        longitudinalIndex(triangleIndex) = longitudinal;

        triangleIndex = triangleIndex+1;
        v1(triangleIndex, :) = vertices(1, :);
        v2(triangleIndex, :) = vertices(3, :);
        v3(triangleIndex, :) = vertices(4, :);
        areas(triangleIndex) = triangleArea(v1(triangleIndex, :), ...
            v2(triangleIndex, :), v3(triangleIndex, :));
        longitudinalIndex(triangleIndex) = longitudinal;
    end
end
assert(all(areas > 0));
end

function area = triangleArea(a, b, c)
area = 0.5*norm(cross(b-a, c-a));
end

function area = sumRoadSurfaceArea(model)
area = sum(model.triangleAreas);
end

function poseBank = buildPoseBank(model, method, param)
candidateCount = numel(model.cumulativeDistance);
poseBank = repmat(struct('groundPoint', zeros(1, 3), ...
    'waypoint', zeros(1, 3), 'opticalAxis', zeros(1, 3), ...
    'groundNormal', zeros(1, 3), 'groundTangent', zeros(1, 3), ...
    'imageAlong', zeros(1, 3), 'imageAcross', zeros(1, 3), ...
    'footprint', zeros(4, 3)), candidateCount, 1);

fixedReference = model.tangent(1, :);
fixedReference(3) = 0;
fixedReference = fixedReference/norm(fixedReference);

for candidateIndex = 1:candidateCount
    ground = model.localCoordinates(candidateIndex, :);
    surfaceNormal = model.normal(candidateIndex, :);
    surfaceTangent = model.tangent(candidateIndex, :);

    if strcmp(method.orientationMode, 'fixed')
        imageAlong = fixedReference- ...
            dot(fixedReference, surfaceNormal)*surfaceNormal;
        assert(norm(imageAlong) > 1e-12, ...
            'Fixed image direction is parallel to the local normal.');
        imageAlong = imageAlong/norm(imageAlong);
    elseif strcmp(method.orientationMode, 'road')
        imageAlong = surfaceTangent- ...
            dot(surfaceTangent, surfaceNormal)*surfaceNormal;
        imageAlong = imageAlong/norm(imageAlong);
    else
        error('Unknown orientation mode: %s', method.orientationMode);
    end
    imageAcross = cross(surfaceNormal, imageAlong);
    imageAcross = imageAcross/norm(imageAcross);

    poseBank(candidateIndex).groundPoint = ground;
    poseBank(candidateIndex).waypoint = ...
        ground+param.cameraDistance*surfaceNormal;
    poseBank(candidateIndex).opticalAxis = -surfaceNormal;
    poseBank(candidateIndex).groundNormal = surfaceNormal;
    poseBank(candidateIndex).groundTangent = surfaceTangent;
    poseBank(candidateIndex).imageAlong = imageAlong;
    poseBank(candidateIndex).imageAcross = imageAcross;
    poseBank(candidateIndex).footprint = tangentPlaneRectangle( ...
        ground, imageAlong, imageAcross, ...
        param.footprintAlong, param.footprintAcross);
end
end

function vertices = tangentPlaneRectangle(centre, along, across, ...
    lengthAlong, lengthAcross)
signs = [-1, -1; 1, -1; 1, 1; -1, 1];
vertices = centre+signs(:, 1)*(lengthAlong/2).*along+ ...
    signs(:, 2)*(lengthAcross/2).*across;
end

function coverageBank = buildCoverageBank(model, poseBank, param)
candidateCount = numel(poseBank);
triangleCount = numel(model.triangleAreas);
areaByTriangle = zeros(candidateCount, triangleCount);
effectiveArea = zeros(candidateCount, 1);

triangleCentres = (model.triangleV1+model.triangleV2+model.triangleV3)/3;
triangleRadius = max([vecnorm(model.triangleV1-triangleCentres, 2, 2), ...
    vecnorm(model.triangleV2-triangleCentres, 2, 2), ...
    vecnorm(model.triangleV3-triangleCentres, 2, 2)], [], 2);
for candidateIndex = 1:candidateCount
    pose = poseBank(candidateIndex);
    offset = triangleCentres-pose.groundPoint;
    alongCoordinate = offset*pose.imageAlong';
    acrossCoordinate = offset*pose.imageAcross';
    possible = abs(alongCoordinate) <= ...
        param.footprintAlong/2+triangleRadius & ...
        abs(acrossCoordinate) <= param.footprintAcross/2+triangleRadius;
    [directions, bounds] = footprintHalfspaces(pose, param);
    for triangleIndex = find(possible(:)')
        areaByTriangle(candidateIndex, triangleIndex) = ...
            clippedTriangleArea(model, triangleIndex, directions, bounds);
    end
    % G(s) is the connected covered road component containing P(s). Road
    % portions that happen to lie under the same image after a distant
    % hairpin are not part of the current observation interval.
    areaByTriangle(candidateIndex, :) = retainAnchorComponent( ...
        areaByTriangle(candidateIndex, :), candidateIndex, model, param);
    effectiveArea(candidateIndex) = sum(areaByTriangle(candidateIndex, :));
end

maximumIndexStep = ceil(param.maximumSpacing/min(model.segmentDistance));
pairOverlap = nan(candidateCount, maximumIndexStep);
pairCoverage = nan(candidateCount, maximumIndexStep);
for firstIndex = 1:candidateCount-1
    for indexStep = 1:min(maximumIndexStep, candidateCount-firstIndex)
        secondIndex = firstIndex+indexStep;
        [overlap, coverage] = exactPairMetrics(firstIndex, secondIndex, ...
            model, poseBank, areaByTriangle, effectiveArea, param);
        pairOverlap(firstIndex, indexStep) = overlap;
        pairCoverage(firstIndex, indexStep) = coverage;
    end
end

coverageBank.areaByTriangle = areaByTriangle;
coverageBank.effectiveArea = effectiveArea;
coverageBank.footprintArea = param.footprintAlong*param.footprintAcross;
coverageBank.maximumIndexStep = maximumIndexStep;
coverageBank.pairOverlap = pairOverlap;
coverageBank.pairCoverage = pairCoverage;
assert(all(effectiveArea > 0));
end

function retainedAreas = retainAnchorComponent( ...
    areas, candidateIndex, model, param)
stripCount = size(model.surfaceX, 1)-1;
stripAreas = accumarray(model.triangleLongitudinalIndex, areas', ...
    [stripCount, 1], @sum, 0);
areaScale = max(model.triangleAreas);
covered = stripAreas > 128*eps(max(1, areaScale));
anchor = min(max(candidateIndex, 1), stripCount);
if ~covered(anchor) && anchor > 1 && covered(anchor-1)
    anchor = anchor-1;
end
assert(covered(anchor), ...
    'The footprint does not cover the road at its anchor point.');
firstStrip = anchor;
while firstStrip > 1 && covered(firstStrip-1)
    firstStrip = firstStrip-1;
end
lastStrip = anchor;
while lastStrip < stripCount && covered(lastStrip+1)
    lastStrip = lastStrip+1;
end
keep = model.triangleLongitudinalIndex >= firstStrip & ...
    model.triangleLongitudinalIndex <= lastStrip;
stripProgress = 0.5*(model.cumulativeDistance(1:end-1)+ ...
    model.cumulativeDistance(2:end));
localWindow = abs(stripProgress-model.cumulativeDistance(candidateIndex)) <= ...
    param.footprintAlong/2+distanceTolerance(model);
keep = keep & localWindow(model.triangleLongitudinalIndex);
retainedAreas = areas;
retainedAreas(~keep') = 0;
end

function [directions, bounds] = footprintHalfspaces(pose, param)
directions = [pose.imageAlong; -pose.imageAlong; ...
    pose.imageAcross; -pose.imageAcross];
bounds = [dot(pose.groundPoint, pose.imageAlong)+param.footprintAlong/2; ...
    dot(pose.groundPoint, -pose.imageAlong)+param.footprintAlong/2; ...
    dot(pose.groundPoint, pose.imageAcross)+param.footprintAcross/2; ...
    dot(pose.groundPoint, -pose.imageAcross)+param.footprintAcross/2];
end

function [overlap, intervalCoverage] = exactPairMetrics( ...
    firstIndex, secondIndex, model, poseBank, areaByTriangle, ...
    effectiveArea, param)
[firstDirections, firstBounds] = footprintHalfspaces( ...
    poseBank(firstIndex), param);
[secondDirections, secondBounds] = footprintHalfspaces( ...
    poseBank(secondIndex), param);
directions = [firstDirections; secondDirections];
bounds = [firstBounds; secondBounds];

firstAreas = areaByTriangle(firstIndex, :);
secondAreas = areaByTriangle(secondIndex, :);
intersectionAreas = zeros(1, numel(model.triangleAreas));
possibleIntersection = firstAreas > 0 & secondAreas > 0;
for triangleIndex = find(possibleIntersection)
    intersectionAreas(triangleIndex) = clippedTriangleArea( ...
        model, triangleIndex, directions, bounds);
end
intersectionArea = sum(intersectionAreas);
overlap = intersectionArea/max(effectiveArea(firstIndex), ...
    effectiveArea(secondIndex));

intervalMask = model.triangleLongitudinalIndex >= firstIndex & ...
    model.triangleLongitudinalIndex < secondIndex;
intervalArea = sum(model.triangleAreas(intervalMask));
unionAreas = firstAreas'+secondAreas'-intersectionAreas';
coveredIntervalArea = sum(unionAreas(intervalMask));
intervalCoverage = min(1, coveredIntervalArea/intervalArea);
end

function area = clippedTriangleArea(model, triangleIndex, directions, bounds)
polygon = [model.triangleV1(triangleIndex, :); ...
    model.triangleV2(triangleIndex, :); ...
    model.triangleV3(triangleIndex, :)];
coordinateScale = max(abs(polygon), [], 'all');
tolerance = 128*eps(max(1, coordinateScale));
for halfspaceIndex = 1:size(directions, 1)
    polygon = clipPolygonByHalfspace(polygon, ...
        directions(halfspaceIndex, :), bounds(halfspaceIndex), tolerance);
    if size(polygon, 1) < 3
        area = 0;
        return;
    end
end
area = 0;
for vertexIndex = 2:size(polygon, 1)-1
    area = area+triangleArea(polygon(1, :), ...
        polygon(vertexIndex, :), polygon(vertexIndex+1, :));
end
end

function output = clipPolygonByHalfspace(polygon, direction, bound, tolerance)
if isempty(polygon)
    output = polygon;
    return;
end
output = zeros(0, 3);
previous = polygon(end, :);
previousValue = dot(previous, direction)-bound;
previousInside = previousValue <= tolerance;
for vertexIndex = 1:size(polygon, 1)
    current = polygon(vertexIndex, :);
    currentValue = dot(current, direction)-bound;
    currentInside = currentValue <= tolerance;
    if currentInside ~= previousInside
        fraction = previousValue/(previousValue-currentValue);
        intersection = previous+fraction*(current-previous);
        output(end+1, :) = intersection; %#ok<AGROW>
    end
    if currentInside
        output(end+1, :) = current; %#ok<AGROW>
    end
    previous = current;
    previousValue = currentValue;
    previousInside = currentInside;
end
end

function [selectedIndices, constraintFeasible] = ...
    selectUniformCandidates(model, coverageBank, param)
candidateCount = numel(model.cumulativeDistance);
largestIndexStep = find(model.cumulativeDistance <= ...
    param.maximumSpacing+distanceTolerance(model), 1, 'last')-1;
assert(largestIndexStep >= 1, ...
    'The candidate interval exceeds the maximum observation spacing.');

for indexStep = largestIndexStep:-1:1
    trial = (1:indexStep:candidateCount)';
    while numel(trial) > 1 && remainingRoadCoverage( ...
            trial(end-1), model, coverageBank) >= ...
            param.minimumRoadCoverage-1e-10
        trial(end) = [];
    end
    if sequenceIsFeasible(trial, model, coverageBank, param) && ...
            remainingRoadCoverage(trial(end), model, coverageBank) >= ...
            param.minimumRoadCoverage-1e-10
        selectedIndices = trial;
        constraintFeasible = true;
        return;
    end
end

% If the fixed-orientation baseline is infeasible even at the finest
% manuscript candidate interval, retain that best-available baseline and
% report the violation explicitly instead of changing any threshold.
selectedIndices = (1:candidateCount)';
while numel(selectedIndices) > 1 && remainingRoadCoverage( ...
        selectedIndices(end-1), model, coverageBank) >= ...
        param.minimumRoadCoverage-1e-10
    selectedIndices(end) = [];
end
overlap = zeros(numel(selectedIndices)-1, 1);
coverage = zeros(numel(selectedIndices)-1, 1);
for intervalIndex = 1:numel(selectedIndices)-1
    [~, overlap(intervalIndex), coverage(intervalIndex)] = ...
        pairIsFeasible(selectedIndices(intervalIndex), ...
        selectedIndices(intervalIndex+1), model, coverageBank, param);
end
[minimumOverlap, overlapIndex] = min(overlap);
[minimumCoverage, coverageIndex] = min(coverage);
constraintFeasible = false;
warning('run_exp1:UniformBaselineInfeasible', ...
    ['No uniform interval satisfies the common constraints. The finest ', ...
    '2 m candidate interval is retained and marked infeasible: minimum ', ...
    'overlap = %.4f at interval %d; minimum interval coverage = %.4f ', ...
    'at interval %d.'], minimumOverlap, overlapIndex, ...
    minimumCoverage, coverageIndex);
end

function selectedIndices = selectAdaptiveCandidates(model, coverageBank, param)
selectedIndices = 1;
currentIndex = 1;
while remainingRoadCoverage(currentIndex, model, coverageBank) < ...
        param.minimumRoadCoverage-1e-10
    candidateIndices = find(model.cumulativeDistance > ...
        model.cumulativeDistance(currentIndex) & ...
        model.cumulativeDistance <= model.cumulativeDistance(currentIndex)+ ...
        param.maximumSpacing+distanceTolerance(model));
    nextIndex = [];
    for candidateIndex = flip(candidateIndices(:)')
        [isFeasible, ~, ~] = pairIsFeasible(currentIndex, ...
            candidateIndex, model, coverageBank, param);
        if isFeasible
            nextIndex = candidateIndex;
            break;
        end
    end
    if isempty(nextIndex)
        error(['No forward candidate from sample %d satisfies the common ', ...
            'overlap and full-coverage constraints.'], currentIndex);
    end
    selectedIndices(end+1, 1) = nextIndex; %#ok<AGROW>
    currentIndex = nextIndex;
end
end

function [localSpacing, farthestIndex] = ...
    computeLocalMaximumFeasibleSpacing(model, coverageBank, param)
% Candidate-wise maximum forward distance under the exact common
% footprint-overlap and interval-coverage constraints. A zero value means
% that no later candidate is feasible from that candidate position.
candidateCount = numel(model.cumulativeDistance);
localSpacing = zeros(candidateCount, 1);
farthestIndex = nan(candidateCount, 1);
for currentIndex = 1:candidateCount-1
    candidateIndices = find(model.cumulativeDistance > ...
        model.cumulativeDistance(currentIndex) & ...
        model.cumulativeDistance <= model.cumulativeDistance(currentIndex)+ ...
        param.maximumSpacing+distanceTolerance(model));
    for candidateIndex = flip(candidateIndices(:)')
        [isFeasible, ~, ~] = pairIsFeasible(currentIndex, ...
            candidateIndex, model, coverageBank, param);
        if isFeasible
            farthestIndex(currentIndex) = candidateIndex;
            localSpacing(currentIndex) = ...
                model.cumulativeDistance(candidateIndex)- ...
                model.cumulativeDistance(currentIndex);
            break;
        end
    end
end
end

function coverage = remainingRoadCoverage(candidateIndex, model, coverageBank)
remainingMask = model.triangleLongitudinalIndex >= candidateIndex;
remainingArea = sum(model.triangleAreas(remainingMask));
if remainingArea <= 0
    coverage = 1;
else
    coveredArea = sum(coverageBank.areaByTriangle( ...
        candidateIndex, remainingMask));
    coverage = min(1, coveredArea/remainingArea);
end
end

function feasible = sequenceIsFeasible(indices, model, coverageBank, param)
feasible = true;
for intervalIndex = 1:numel(indices)-1
    [pairFeasible, ~, ~] = pairIsFeasible(indices(intervalIndex), ...
        indices(intervalIndex+1), model, coverageBank, param);
    if ~pairFeasible
        feasible = false;
        return;
    end
end
end

function [feasible, overlap, intervalCoverage] = pairIsFeasible( ...
    firstIndex, secondIndex, model, coverageBank, param)
assert(secondIndex > firstIndex);
distance = model.cumulativeDistance(secondIndex)- ...
    model.cumulativeDistance(firstIndex);
indexStep = secondIndex-firstIndex;
assert(indexStep <= coverageBank.maximumIndexStep);
overlap = coverageBank.pairOverlap(firstIndex, indexStep);
intervalCoverage = coverageBank.pairCoverage(firstIndex, indexStep);

numericTolerance = 1e-10;
feasible = distance <= param.maximumSpacing+distanceTolerance(model) && ...
    overlap+numericTolerance >= param.minimumOverlap && ...
    intervalCoverage+numericTolerance >= param.minimumRoadCoverage;
end

function tolerance = distanceTolerance(model)
tolerance = 128*eps(max(1, model.cumulativeDistance(end)));
end

function result = assembleResult(model, poseBank, coverageBank, ...
    selectedIndices, param)
poses = poseBank(selectedIndices);
count = numel(selectedIndices);
groundPoints = reshape([poses.groundPoint], 3, count)';
waypoints = reshape([poses.waypoint], 3, count)';
opticalAxes = reshape([poses.opticalAxis], 3, count)';
groundNormals = reshape([poses.groundNormal], 3, count)';
groundTangents = reshape([poses.groundTangent], 3, count)';
imageAlong = reshape([poses.imageAlong], 3, count)';
imageAcross = reshape([poses.imageAcross], 3, count)';
footprints = repmat(struct('vertices', zeros(4, 3)), count, 1);
for waypointIndex = 1:count
    footprints(waypointIndex).vertices = poses(waypointIndex).footprint;
end

adjacentOverlap = zeros(max(0, count-1), 1);
intervalCoverage = zeros(max(0, count-1), 1);
for intervalIndex = 1:count-1
    [~, adjacentOverlap(intervalIndex), intervalCoverage(intervalIndex)] = ...
        pairIsFeasible(selectedIndices(intervalIndex), ...
        selectedIndices(intervalIndex+1), model, coverageBank, param);
end

coveredRoadArea = 0;
for intervalIndex = 1:count-1
    intervalMask = model.triangleLongitudinalIndex >= ...
        selectedIndices(intervalIndex) & ...
        model.triangleLongitudinalIndex < selectedIndices(intervalIndex+1);
    intervalArea = sum(model.triangleAreas(intervalMask));
    coveredRoadArea = coveredRoadArea+ ...
        intervalCoverage(intervalIndex)*intervalArea;
end
remainingMask = model.triangleLongitudinalIndex >= selectedIndices(end);
coveredRoadArea = coveredRoadArea+sum(coverageBank.areaByTriangle( ...
    selectedIndices(end), remainingMask));
totalRoadCoverage = min(1, coveredRoadArea/sum(model.triangleAreas));
footprintUtilization = coverageBank.effectiveArea(selectedIndices)/ ...
    coverageBank.footprintArea;

result.selectedIndices = selectedIndices;
result.routeProgress = model.cumulativeDistance(selectedIndices);
result.groundPoints = groundPoints;
result.waypoints = waypoints;
result.opticalAxes = opticalAxes;
result.groundNormals = groundNormals;
result.groundTangents = groundTangents;
result.imageAlong = imageAlong;
result.imageAcross = imageAcross;
result.footprints = footprints;
result.intervalDistance = diff(result.routeProgress);
result.adjacentOverlap = adjacentOverlap;
result.intervalCoverage = intervalCoverage;
result.totalRoadCoverage = totalRoadCoverage;
result.footprintUtilization = footprintUtilization;
end

function validateResult(result, model, method, param)
assert(all(isfinite(result.waypoints), 'all'));
assert(all(isfinite(result.opticalAxes), 'all'));
assert(result.selectedIndices(1) == 1);
assert(all(diff(result.selectedIndices) > 0));
assert(all(result.selectedIndices == round(result.selectedIndices)));
assert(all(result.intervalDistance <= ...
    param.maximumSpacing+distanceTolerance(model)));
if result.constraintFeasible
    assert(all(result.adjacentOverlap+1e-10 >= param.minimumOverlap));
end
assert(all(result.intervalCoverage+1e-12 >= ...
    param.minimumRoadCoverage));
assert(result.totalRoadCoverage+1e-12 >= param.minimumRoadCoverage);
if strcmp(method.samplingMode, 'adaptive')
    assert(result.constraintFeasible);
    selectedStarts = result.selectedIndices(1:end-1);
    expectedSpacing = result.localMaximumFeasibleSpacing(selectedStarts);
    assert(all(expectedSpacing > 0));
    assert(max(abs(expectedSpacing-result.intervalDistance)) < ...
        distanceTolerance(model));
end

cameraDistance = vecnorm(result.waypoints-result.groundPoints, 2, 2);
assert(max(abs(cameraDistance-param.cameraDistance)) < 1e-9);
assert(max(vecnorm(result.opticalAxes+result.groundNormals, 2, 2)) ...
    < 1e-10);
assert(max(abs(sum(result.imageAlong.*result.groundNormals, 2))) < 1e-10);
assert(max(abs(sum(result.imageAcross.*result.groundNormals, 2))) < 1e-10);

for footprintIndex = 1:numel(result.footprints)
    vertices = result.footprints(footprintIndex).vertices;
    sideLengths = vecnorm(vertices([2, 3, 4, 1], :)-vertices, 2, 2);
    assert(abs(sideLengths(1)-param.footprintAlong) < 1e-8);
    assert(abs(sideLengths(2)-param.footprintAcross) < 1e-8);
end

if strcmp(method.orientationMode, 'road')
    alignment = abs(sum(result.imageAlong.*result.groundTangents, 2));
    assert(min(alignment) > 1-1e-10);
elseif strcmp(method.orientationMode, 'fixed')
    fixedReference = model.tangent(1, :);
    fixedReference(3) = 0;
    fixedReference = fixedReference/norm(fixedReference);
    expected = fixedReference- ...
        sum(result.groundNormals.*fixedReference, 2).*result.groundNormals;
    expected = expected./vecnorm(expected, 2, 2);
    assert(max(vecnorm(result.imageAlong-expected, 2, 2)) < 1e-10);
end

if strcmp(method.samplingMode, 'uniform') && numel(result.selectedIndices) > 2
    indexSteps = diff(result.selectedIndices);
    assert(isscalar(unique(indexSteps(1:end-1))));
end
end

function comparison = findComparisonLocation(model, roadResults, methods, param)
fixedReference = model.tangent(1, :);
fixedReference(3) = 0;
fixedReference = fixedReference/norm(fixedReference);
fixedProjection = fixedReference- ...
    sum(model.normal.*fixedReference, 2).*model.normal;
fixedProjection = fixedProjection./vecnorm(fixedProjection, 2, 2);
directionMismatch = acos(max(-1, min(1, ...
    abs(sum(fixedProjection.*model.tangent, 2)))));

edgeGuardDistance = param.footprintAlong/2;
interior = model.cumulativeDistance >= edgeGuardDistance & ...
    model.cumulativeDistance <= model.cumulativeDistance(end)-edgeGuardDistance;
directionMismatch(~interior) = -inf;
[maximumMismatch, sampleIndex] = max(directionMismatch);

comparison.sampleIndex = sampleIndex;
comparison.routeProgress = model.cumulativeDistance(sampleIndex);
comparison.directionMismatchDeg = rad2deg(maximumMismatch);
comparison.windowHalfLength = param.footprintAlong;
comparison.roadWidthScale = 6; % Display only; computations use true width.
comparison.nearestWaypoint = zeros(1, numel(methods));
for methodIndex = 1:numel(methods)
    [~, comparison.nearestWaypoint(methodIndex)] = min(abs( ...
        roadResults(methodIndex).routeProgress-comparison.routeProgress));
end
end

function drawRoadComparisonFigure(model, roadResults, methods, ...
    comparison, outputDir)
style.figureSize = [40, 80, 2400, 1000];
style.roadEdge = [0.08, 0.08, 0.08];
style.centreline = [0.04, 0.04, 0.04];
style.fontSize = 15;
style.panelFontSize = 17;
style.metricFontSize = 14;
style.footprintAlpha = 0.045;
style.footprintEdgeAlpha = 0.34;
style.arrowLength = 7.5;
style.verticalExaggeration = 3.3;

[displayX, displayY, displayZ, displayElevation] = ...
    expandedRoadSurface(model, comparison.roadWidthScale);
[commonLimits, commonElevationLimits] = commonFigureLimits( ...
    displayX, displayY, displayZ, displayElevation, roadResults);

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', style.figureSize);
set(fig, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman');
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', ...
    'Padding', 'loose');
axesHandles = gobjects(1, 3);

for methodIndex = 1:numel(methods)
    ax = nexttile(layout, methodIndex);
    axesHandles(methodIndex) = ax;
    drawRoadPanel(ax, model, roadResults(methodIndex), ...
        methods(methodIndex), displayX, displayY, displayZ, ...
        displayElevation, commonLimits, commonElevationLimits, style);
end
drawnow;

for methodIndex = 1:numel(methods)
    position = axesHandles(methodIndex).Position;
    panelLabel = sprintf('(%s) Road %s: %s', methods(methodIndex).id, ...
        model.routeId, methods(methodIndex).name);
    annotation(fig, 'textbox', ...
        [position(1)+0.006, position(2)+position(4)-0.043, 0.23, 0.036], ...
        'String', panelLabel, 'FitBoxToText', 'on', 'EdgeColor', 'none', ...
        'BackgroundColor', 'w', 'Margin', 1.5, ...
        'FontName', 'Times New Roman', 'FontSize', style.panelFontSize, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    result = roadResults(methodIndex);
    metricLabel = sprintf(['N = %d\nCoverage = %.1f%%\n', ...
        'Min. overlap = %.1f%%'], numel(result.selectedIndices), ...
        100*result.totalRoadCoverage, 100*min(result.adjacentOverlap));
    annotation(fig, 'textbox', ...
        [position(1)+0.012, position(2)+0.050, 0.13, 0.085], ...
        'String', metricLabel, 'FitBoxToText', 'on', ...
        'EdgeColor', [0.30, 0.30, 0.30], 'BackgroundColor', 'w', ...
        'FaceAlpha', 0.90, 'Margin', 4, ...
        'FontName', 'Times New Roman', 'FontSize', style.metricFontSize, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

    drawComparisonInset(fig, axesHandles(methodIndex), model, ...
        result, methods(methodIndex), methodIndex, comparison, style);
end

colorbarObject = colorbar(axesHandles(3), 'eastoutside');
colorbarObject.Layout.Tile = 'east';
colorbarObject.Label.String = 'Elevation (m)';
colorbarObject.Label.FontName = 'Times New Roman';
colorbarObject.Label.FontSize = style.fontSize;
colorbarObject.FontName = 'Times New Roman';
colorbarObject.FontSize = style.fontSize;
drawnow;

fileStem = sprintf('fig_exp1_road_%s_comparison', model.routeId);
pngPath = fullfile(outputDir, [fileStem, '.png']);
pdfPath = fullfile(outputDir, [fileStem, '.pdf']);
exportgraphics(fig, pngPath, 'Resolution', 300, ...
    'BackgroundColor', 'white');
exportgraphics(fig, pdfPath, 'ContentType', 'image', ...
    'Resolution', 600, 'BackgroundColor', 'white');
close(fig);
end

function drawRoadPanel(ax, model, result, method, ...
    displayX, displayY, displayZ, displayElevation, ...
    commonLimits, commonElevationLimits, style)
hold(ax, 'on');
surf(ax, displayX, displayY, displayZ, displayElevation, ...
    'FaceColor', 'interp', 'FaceAlpha', 0.98, 'EdgeColor', 'none');
plot3(ax, displayX(:, 1), displayY(:, 1), displayZ(:, 1), '-', ...
    'Color', style.roadEdge, 'LineWidth', 0.9, 'HandleVisibility', 'off');
plot3(ax, displayX(:, end), displayY(:, end), displayZ(:, end), '-', ...
    'Color', style.roadEdge, 'LineWidth', 0.9, 'HandleVisibility', 'off');
plot3(ax, model.localCoordinates(:, 1), model.localCoordinates(:, 2), ...
    model.localCoordinates(:, 3), '-', 'Color', style.centreline, ...
    'LineWidth', 1.25, 'HandleVisibility', 'off');

for footprintIndex = 1:numel(result.footprints)
    vertices = result.footprints(footprintIndex).vertices+ ...
        0.05*result.groundNormals(footprintIndex, :);
    patch(ax, 'Vertices', vertices, 'Faces', 1:4, ...
        'FaceColor', method.color, 'FaceAlpha', style.footprintAlpha, ...
        'EdgeColor', method.color, 'EdgeAlpha', style.footprintEdgeAlpha, ...
        'LineWidth', 0.52, 'HandleVisibility', 'off');
end

quiver3(ax, result.waypoints(:, 1), result.waypoints(:, 2), ...
    result.waypoints(:, 3), style.arrowLength*result.opticalAxes(:, 1), ...
    style.arrowLength*result.opticalAxes(:, 2), ...
    style.arrowLength*result.opticalAxes(:, 3), 0, ...
    'Color', method.color, 'LineWidth', 0.78, 'MaxHeadSize', 0.60, ...
    'HandleVisibility', 'off');
plot3(ax, result.waypoints(:, 1), result.waypoints(:, 2), ...
    result.waypoints(:, 3), 'o', 'Color', method.color, ...
    'MarkerFaceColor', method.color, 'MarkerSize', 3.8, ...
    'LineWidth', 0.4, 'HandleVisibility', 'off');

xlabel(ax, 'x (m)');
ylabel(ax, 'y (m)');
zlabel(ax, '\Deltaz (m)');
set(ax, 'FontName', 'Times New Roman', 'FontSize', style.fontSize, ...
    'LineWidth', 0.8, 'Box', 'on', 'Layer', 'top', ...
    'GridColor', [0.72, 0.72, 0.72], 'GridAlpha', 0.30);
grid(ax, 'on');
xlim(ax, commonLimits(1, :));
ylim(ax, commonLimits(2, :));
zlim(ax, commonLimits(3, :));
clim(ax, commonElevationLimits);
daspect(ax, [1, 1, 1/style.verticalExaggeration]);
view(ax, 36, 27);
camproj(ax, 'orthographic');
colormap(ax, parula(256));
end

function drawComparisonInset(fig, mainAxes, model, result, method, ...
    methodIndex, comparison, style)
mainPosition = mainAxes.Position;
insetPosition = [mainPosition(1)+0.58*mainPosition(3), ...
    mainPosition(2)+0.57*mainPosition(4), ...
    0.37*mainPosition(3), 0.36*mainPosition(4)];
insetAxes = axes(fig, 'Position', insetPosition);
hold(insetAxes, 'on');

[displayX, displayY, displayZ, displayElevation] = ...
    expandedRoadSurface(model, comparison.roadWidthScale);
roadMask = abs(model.cumulativeDistance-comparison.routeProgress) <= ...
    comparison.windowHalfLength;
surf(insetAxes, displayX(roadMask, :), displayY(roadMask, :), ...
    displayZ(roadMask, :), displayElevation(roadMask, :), ...
    'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.98);
plot3(insetAxes, displayX(roadMask, 1), displayY(roadMask, 1), ...
    displayZ(roadMask, 1), '-', 'Color', style.roadEdge, 'LineWidth', 0.8);
plot3(insetAxes, displayX(roadMask, end), displayY(roadMask, end), ...
    displayZ(roadMask, end), '-', 'Color', style.roadEdge, 'LineWidth', 0.8);
plot3(insetAxes, model.localCoordinates(roadMask, 1), ...
    model.localCoordinates(roadMask, 2), model.localCoordinates(roadMask, 3), ...
    '-', 'Color', style.centreline, 'LineWidth', 1.1);

nearFootprints = abs(result.routeProgress-comparison.routeProgress) <= ...
    comparison.windowHalfLength/2;
for footprintIndex = find(nearFootprints(:)')
    vertices = result.footprints(footprintIndex).vertices+ ...
        0.05*result.groundNormals(footprintIndex, :);
    patch(insetAxes, 'Vertices', vertices, 'Faces', 1:4, ...
        'FaceColor', method.color, 'FaceAlpha', 0.10, ...
        'EdgeColor', method.color, 'EdgeAlpha', 0.65, 'LineWidth', 0.8);
end

waypointIndex = comparison.nearestWaypoint(methodIndex);
arrowLength = style.arrowLength;
quiver3(insetAxes, result.waypoints(waypointIndex, 1), ...
    result.waypoints(waypointIndex, 2), result.waypoints(waypointIndex, 3), ...
    arrowLength*result.opticalAxes(waypointIndex, 1), ...
    arrowLength*result.opticalAxes(waypointIndex, 2), ...
    arrowLength*result.opticalAxes(waypointIndex, 3), 0, ...
    'Color', method.color, 'LineWidth', 1.3, 'MaxHeadSize', 0.55);
plot3(insetAxes, result.waypoints(nearFootprints, 1), ...
    result.waypoints(nearFootprints, 2), result.waypoints(nearFootprints, 3), ...
    'o', 'Color', method.color, 'MarkerFaceColor', method.color, ...
    'MarkerSize', 4.2);

insetPoints = [reshape(displayX(roadMask, :), [], 1), ...
    reshape(displayY(roadMask, :), [], 1), ...
    reshape(displayZ(roadMask, :), [], 1)];
for footprintIndex = find(nearFootprints(:)')
    insetPoints = [insetPoints; ...
        result.footprints(footprintIndex).vertices]; %#ok<AGROW>
end
insetPoints = [insetPoints; result.waypoints(nearFootprints, :)];
minimum = min(insetPoints, [], 1);
maximum = max(insetPoints, [], 1);
margin = max(0.06*(maximum-minimum), [1.5, 1.5, 1.5]);
xlim(insetAxes, [minimum(1)-margin(1), maximum(1)+margin(1)]);
ylim(insetAxes, [minimum(2)-margin(2), maximum(2)+margin(2)]);
zlim(insetAxes, [minimum(3)-margin(3), maximum(3)+margin(3)]);
clim(insetAxes, [min(displayElevation(roadMask, :), [], 'all'), ...
    max(displayElevation(roadMask, :), [], 'all')]);
view(insetAxes, 36, 27);
camproj(insetAxes, 'orthographic');
daspect(insetAxes, [1, 1, 1/style.verticalExaggeration]);
colormap(insetAxes, parula(256));
set(insetAxes, 'XTick', [], 'YTick', [], 'ZTick', [], ...
    'Box', 'on', 'LineWidth', 1.0, 'Color', 'w');
end

function [limits, elevationLimits] = commonFigureLimits( ...
    surfaceX, surfaceY, surfaceZ, surfaceElevation, roadResults)
points = [surfaceX(:), surfaceY(:), surfaceZ(:)];
for methodIndex = 1:numel(roadResults)
    points = [points; roadResults(methodIndex).waypoints]; %#ok<AGROW>
    for footprintIndex = 1:numel(roadResults(methodIndex).footprints)
        points = [points; ...
            roadResults(methodIndex).footprints(footprintIndex).vertices]; ...
            %#ok<AGROW>
    end
end
minimum = min(points, [], 1);
maximum = max(points, [], 1);
margin = max(0.035*(maximum-minimum), [2, 2, 2]);
limits = [minimum(1)-margin(1), maximum(1)+margin(1); ...
    minimum(2)-margin(2), maximum(2)+margin(2); ...
    minimum(3)-margin(3), maximum(3)+margin(3)];
elevationLimits = [min(surfaceElevation, [], 'all'), ...
    max(surfaceElevation, [], 'all')];
if diff(elevationLimits) < eps(max(abs(elevationLimits)))
    elevationLimits = elevationLimits+[-0.5, 0.5];
end
end

function [surfaceX, surfaceY, surfaceZ, surfaceElevation] = ...
    expandedRoadSurface(model, widthScale)
offsets = linspace(-widthScale*model.width/2, ...
    widthScale*model.width/2, 9);
surfaceX = model.localCoordinates(:, 1)+model.lateral(:, 1)*offsets;
surfaceY = model.localCoordinates(:, 2)+model.lateral(:, 2)*offsets;
surfaceZ = model.localCoordinates(:, 3)+model.lateral(:, 3)*offsets;
surfaceElevation = model.worldCoordinates(:, 3)+model.lateral(:, 3)*offsets;
end
