function run_exp2b()
%RUN_EXP2B Plot representative GA-DP trajectories for R = 5, 10, and 20.
% Run this file directly in MATLAB. It reads the measured-road waypoint
% library produced by run_exp2a, selects turn-rich real-road tasks, solves each
% task with the same GA-DP settings, and exports a three-panel manuscript
% figure. No synthetic roads, elevations, or waypoints are generated.

close all;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
projectRoot = fileparts(fileparts(scriptDir));
dataDir = fullfile(projectRoot, 'dem');
outputDir = fullfile(projectRoot, 'results', 'exp2');
inputFile = fullfile(outputDir, 'exp2a_results.mat');
if ~isfolder(outputDir)
    mkdir(outputDir);
end
if ~isfile(inputFile)
    addpath(scriptDir);
    run_exp2a();
end
assert(isfile(inputFile), 'Experiment-2a results not found: %s', inputFile);
storedVariables = string({whos('-file',inputFile).name});
if ~ismember("masterProblem",storedVariables)
    run_exp2a();
end
assert(exist('shaperead', 'file') == 2 && ...
    exist('readgeoraster', 'file') == 2, ...
    'Mapping Toolbox is required to draw the measured satellite map.');

input = load(inputFile, 'routeLibrary', 'instances', 'launchPoint', ...
    'recoveryPoint', 'param', 'cfg','masterProblem');
routeLibrary = input.routeLibrary;
launchPoint = input.launchPoint;
recoveryPoint = input.recoveryPoint;
param = input.param;
ga = input.cfg.ga;

%% Figure experiment settings
cfg.scales = [5, 10, 20];
cfg.independentRuns = 20;
cfg.seed = 20264000;
% Per-scale instances with the greatest certified lateral deviation and
% heading-change count in the exact solutions. Quantitative results still
% aggregate all ten instances at every scale.
cfg.instanceNumberByScale = [5,2,2];
cfg.outputStem = 'fig_exp2b_nested_gadp_trajectories';
cfg.selectedColor = [0.00, 0.31, 0.76];
cfg.transferColor = [0.95, 0.48, 0.08];
cfg.launchColor = [0.10, 0.62, 0.24];
cfg.recoveryColor = [0.88, 0.13, 0.10];
cfg.directionColor = [1.00, 0.82, 0.00];
cfg.detourColor = [0.90, 0.05, 0.55];
cfg.unsafeChordColor = [0.72, 0.72, 0.72];
cfg.visibleDetourThreshold_m = 20;
cfg.targetLineWidth = 1.35;
cfg.targetHaloWidth = 2.25;
cfg.transferLineWidth = 0.70;
cfg.transferHaloWidth = 1.35;
cfg.fontSize = 7.5;
cfg.panelFontSize = 7.6;
cfg.metricFontSize = 6.8;
cfg.showMetrics = false;
cfg.routeLabelFontSize = 6.2;
cfg.figureInches = [0.5, 0.5, 7.5, 3.05];
cfg.mapMarginFraction = 0.12;

%% Deterministic scale-specific representative tasks
tasks = repmat(struct('R', 0, 'instanceNumber',0,'libraryIndices', [], ...
    'solution', [], 'problem', [], 'bestRun', 0, 'bestSeed', 0), 1, 3);
for taskIndex = 1:numel(cfg.scales)
    scale = cfg.scales(taskIndex);
    instanceNumber = cfg.instanceNumberByScale(taskIndex);
    match = find([input.instances.R] == scale & ...
        [input.instances.number] == instanceNumber);
    assert(isscalar(match));
    tasks(taskIndex).R = scale;
    tasks(taskIndex).instanceNumber = instanceNumber;
    tasks(taskIndex).libraryIndices = input.instances(match).libraryIndices;
    assert(numel(tasks(taskIndex).libraryIndices) == scale && ...
        ~ismember(1,tasks(taskIndex).libraryIndices), ...
        'The boundary corridor must be excluded from target roads.');
end

%% Solve every representative task with GA-DP and retain its best trajectory
summaryRows = repmat(struct('ScaleR', 0, 'RoadCount', 0, ...
    'IndependentRuns', 0, 'BestRun', 0, 'BestSeed', 0, ...
    'J_E_kJ', NaN, 'L3D_m', NaN, 'Hplus_m', NaN, ...
    'FlightTime_s',NaN,'TransferEnergy_kJ',NaN, ...
    'DetourCount',NaN,'UnsafeStraightCount',NaN, ...
    'OrderLibraryIndices', "", 'OrderRouteIDs', "", ...
    'Directions', ""), 0, 1);
for taskIndex = 1:numel(tasks)
    problem = subsetRoutingProblem(input.masterProblem, ...
        tasks(taskIndex).libraryIndices);
    heuristicOrder = nearestNeighbourOrder(problem);
    heuristicSolution = solveFixedOrderDP(heuristicOrder, problem);
    bestSolution = [];
    bestRun = 0;
    bestSeed = 0;
    fprintf('Representative R=%d instance %d: %d independent GA-DP runs\n', ...
        tasks(taskIndex).R,tasks(taskIndex).instanceNumber,cfg.independentRuns);
    for runIndex = 1:cfg.independentRuns
        seed = cfg.seed+1000*tasks(taskIndex).R+runIndex;
        solution = solveGaDp(problem, ga, seed, heuristicSolution);
        validateSolution(solution, problem);
        if isempty(bestSolution) || solution.J < bestSolution.J
            bestSolution = solution;
            bestRun = runIndex;
            bestSeed = seed;
        end
    end
    tasks(taskIndex).solution = bestSolution;
    tasks(taskIndex).problem = problem;
    tasks(taskIndex).bestRun = bestRun;
    tasks(taskIndex).bestSeed = bestSeed;

    globalOrder = tasks(taskIndex).libraryIndices(bestSolution.order);
    symbols = repmat("+", 1, tasks(taskIndex).R);
    symbols(bestSolution.directions == 2) = "-";
    row.ScaleR = tasks(taskIndex).R;
    row.RoadCount = tasks(taskIndex).R;
    row.IndependentRuns = cfg.independentRuns;
    row.BestRun = bestRun;
    row.BestSeed = bestSeed;
    row.J_E_kJ = bestSolution.J;
    row.L3D_m = bestSolution.L3D;
    row.Hplus_m = bestSolution.Hplus;
    row.FlightTime_s = bestSolution.Time_s;
    row.TransferEnergy_kJ = bestSolution.TransferEnergy_kJ;
    row.DetourCount = bestSolution.DetourCount;
    row.UnsafeStraightCount = bestSolution.UnsafeDirectCount;
    row.OrderLibraryIndices = join(string(globalOrder), '-');
    row.OrderRouteIDs = join(string({routeLibrary(globalOrder).routeId}), '|');
    row.Directions = join(symbols, '');
    summaryRows(end+1, 1) = row; %#ok<AGROW>
    fprintf('  best run %d: J=%.2f, L3D=%.2f m, Hplus=%.2f m\n', ...
        bestRun, bestSolution.J, bestSolution.L3D, bestSolution.Hplus);
end
summary = struct2table(summaryRows);

%% Reconstruct the measured ground centrelines used by the map figure
majorPath = findInput(dataDir, {'major_roads.shp', '大路.shp'});
minorPath = findInput(dataDir, {'minor_roads.shp', '小路.shp'});
imagePath = findInput(dataDir, {'orthophoto.tif', 'cgs2000.tif'});
majorShapes = shaperead(majorPath);
minorShapes = shaperead(minorPath);
centrelineByLibrary = cell(1, numel(routeLibrary));
figureIndices = unique([tasks.libraryIndices]);
for libraryIndex = figureIndices
    centrelineByLibrary{libraryIndex} = reconstructCentreline( ...
        routeLibrary(libraryIndex), majorShapes, minorShapes, param);
end

%% Common satellite crop for direct visual comparison across scales
[orthophoto, imageReference] = readgeoraster(imagePath);
if size(orthophoto, 3) > 3
    orthophoto = orthophoto(:, :, 1:3);
end
[xLimits, yLimits] = commonMapLimits(figureIndices, ...
    centrelineByLibrary,launchPoint,recoveryPoint,imageReference, ...
    cfg.mapMarginFraction);
[imageCrop, cropReference] = cropMap( ...
    orthophoto, imageReference, xLimits, yLimits);
clear orthophoto;

%% Three-panel cross-column figure
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', cfg.figureInches);
set(fig, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman');
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
axesHandles = gobjects(1, 3);
legendHandles = gobjects(1, 6);

for panelIndex = 1:3
    ax = nexttile(layout, panelIndex);
    axesHandles(panelIndex) = ax;
    mapshow(imageCrop, cropReference, 'Parent', ax);
    hold(ax, 'on');
    [targetHandle, transferHandle, detourHandle, launchHandle, ...
        recoveryHandle, directionHandle] = ...
        drawTrajectoryPanel(ax, tasks(panelIndex), routeLibrary, ...
        centrelineByLibrary,launchPoint,recoveryPoint,cfg);
    if panelIndex == 1
        legendHandles = [targetHandle, transferHandle,detourHandle, ...
            launchHandle, recoveryHandle, directionHandle];
    end

    xlim(ax, xLimits);
    ylim(ax, yLimits);
    axis(ax, 'equal');
    set(ax, 'YDir', 'normal', 'FontName', 'Times New Roman', ...
        'FontSize', cfg.fontSize, 'LineWidth', 0.60, ...
        'Box', 'on', 'Layer', 'top', 'TickDir', 'out');
    ax.XAxis.Exponent = 5;
    ax.YAxis.Exponent = 6;
    xlabel(ax, 'Easting (m)');
    if panelIndex == 1
        ylabel(ax, 'Northing (m)');
    else
        yticklabels(ax, []);
    end
    panelText = sprintf('(%c)  R = %d', 'a'+panelIndex-1, ...
        tasks(panelIndex).R);
    text(ax, 0.023, 0.977, panelText, 'Units', 'normalized', ...
        'Color', [0.08, 0.08, 0.08], 'FontWeight', 'bold', ...
        'FontSize', cfg.panelFontSize, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'Clipping', 'on');
    text(ax, 0.020, 0.980, panelText, 'Units', 'normalized', ...
        'Color', 'w', 'FontWeight', 'bold', ...
        'FontSize', cfg.panelFontSize, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'Clipping', 'on');

    if cfg.showMetrics
        solution = tasks(panelIndex).solution;
        metricText = {sprintf('$J=%.1f\\,\\mathrm{kJ}$', solution.J), ...
            sprintf('$L_{3D}=%.1f\\,\\mathrm{m}$', solution.L3D), ...
            sprintf('$H^{+}=%.1f\\,\\mathrm{m}$', solution.Hplus)};
        text(ax, 0.978, 0.032, metricText, 'Units', 'normalized', ...
            'Interpreter', 'latex', 'FontSize', cfg.metricFontSize, ...
            'Color', [0.04, 0.04, 0.04], ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        text(ax, 0.975, 0.035, metricText, 'Units', 'normalized', ...
            'Interpreter', 'latex', 'FontSize', cfg.metricFontSize, ...
            'Color', 'w', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end
end

legendObject = legend(axesHandles(1), legendHandles, ...
    {'Target road', 'Inter-road flight','Energy-optimal graph path', 'Launch point', ...
        'Recovery point', 'Visit direction'}, ...
    'Orientation', 'horizontal', 'NumColumns', 6, ...
    'Box', 'off', 'FontName', 'Times New Roman', ...
    'FontSize', cfg.fontSize, 'Location', 'southoutside');
legendObject.Layout.Tile = 'south';

pngPath = fullfile(outputDir, [cfg.outputStem, '.png']);
pdfPath = fullfile(outputDir, [cfg.outputStem, '.pdf']);
exportgraphics(fig, pngPath, 'Resolution', 400, ...
    'BackgroundColor', 'white');
exportgraphics(fig, pdfPath, 'ContentType', 'image', ...
    'Resolution', 400, 'BackgroundColor', 'white');
close(fig);

writetable(summary, fullfile(outputDir, 'exp2b_summary.csv'), ...
    'Encoding', 'UTF-8');
metadata.generatedAt = datetime('now', 'TimeZone', 'Asia/Shanghai');
metadata.sourceExperiment = 'exp2a_results.mat';
metadata.representativeTasks = string(arrayfun(@(scale,instance) ...
    sprintf('R%02d_I%02d',scale,instance),cfg.scales, ...
    cfg.instanceNumberByScale,'UniformOutput',false));
metadata.selectionRule = ['maximum visible lateral-deviation and ', ...
    'heading-change evidence within each scale'];
metadata.boundaryDefinition = ['Distinct fixed launch and recovery points ', ...
    'at opposite ends of a measured boundary corridor excluded from all tasks'];
save(fullfile(outputDir, 'exp2b_results.mat'), 'tasks', 'summary', ...
    'centrelineByLibrary', 'launchPoint', 'recoveryPoint', 'cfg', 'metadata', ...
    'xLimits', 'yLimits');

finalAudit(tasks,summary,routeLibrary,launchPoint,recoveryPoint,xLimits,yLimits, ...
    pngPath, pdfPath);
fprintf('Experiment 2b completed. Figure written to:\n  %s\n  %s\n', ...
    pngPath, pdfPath);
end

%% Legacy spatial-subset helper retained for reproducibility

function subset = selectSpatiallyNestedSubset( ...
    parentIndices, targetCount, routeLibrary, launchPoint)
assert(targetCount <= numel(parentIndices));
coordinates = zeros(numel(parentIndices), 2);
for index = 1:numel(parentIndices)
    waypoints = routeLibrary(parentIndices(index)).waypoints;
    coordinates(index, :) = mean(waypoints(:, 1:2), 1);
end
[~,anchor] = min(vecnorm(coordinates-launchPoint(1:2),2,2));
selected = false(1, numel(parentIndices));
selected(anchor) = true;
minimumDistance = vecnorm(coordinates-coordinates(anchor, :), 2, 2);
while sum(selected) < targetCount
    score = minimumDistance;
    score(selected) = -inf;
    [~, next] = max(score);
    selected(next) = true;
    distance = vecnorm(coordinates-coordinates(next, :), 2, 2);
    minimumDistance = min(minimumDistance, distance);
end
subset = parentIndices(selected);
end

%% Routing cost and GA-DP

function problem = subsetRoutingProblem(master,indices)
problem.R = numel(indices);
problem.launchPoint = master.launchPoint;
problem.recoveryPoint = master.recoveryPoint;
two = 1:2;
fields2 = {'entry','exit','entryVelocity','exitVelocity'};
for k=1:numel(fields2)
    name=fields2{k}; problem.(name)=master.(name)(indices,two,:);
end
fields = {'internalL','internalH','internalJ','internalTime', ...
    'depotEntryL','depotEntryH','depotEntryJ','depotEntryTime', ...
    'depotEntryPath','depotEntryMinimumClearance', ...
    'depotEntryDirectMinimumClearance','depotEntryDetour', ...
    'exitDepotL','exitDepotH','exitDepotJ','exitDepotTime', ...
    'exitDepotPath','exitDepotMinimumClearance', ...
    'exitDepotDirectMinimumClearance','exitDepotDetour'};
for k=1:numel(fields)
    name=fields{k}; problem.(name)=master.(name)(indices,two);
end
fields4 = {'transitionL','transitionH','transitionJ','transitionTime', ...
    'transitionPath','transitionMinimumClearance', ...
    'transitionDirectMinimumClearance','transitionDetour'};
for k=1:numel(fields4)
    name=fields4{k};
    problem.(name)=master.(name)(indices,two,indices,two);
end
problem.requiredClearance=master.requiredClearance;
end

function problem = buildRoutingProblem(routes,launchPoint,recoveryPoint,param)
R = numel(routes);
entry = zeros(R, 2, 3);
exitPoint = zeros(R, 2, 3);
internalL = zeros(R, 2);
internalH = zeros(R, 2);
internalJ = zeros(R, 2);
for road = 1:R
    waypoints = routes(road).waypoints;
    for direction = 1:2
        if direction == 1
            path = waypoints;
        else
            path = flipud(waypoints);
        end
        entry(road, direction, :) = path(1, :);
        exitPoint(road, direction, :) = path(end, :);
        for waypoint=1:size(path,1)-1
            [energy_kJ,length3D,climb] = legMetrics( ...
                path(waypoint,:),path(waypoint+1,:),param);
            internalJ(road,direction)=internalJ(road,direction)+energy_kJ;
            internalL(road,direction)=internalL(road,direction)+length3D;
            internalH(road,direction)=internalH(road,direction)+climb;
        end
    end
end

problem.R = R;
problem.internalL = internalL;
problem.internalH = internalH;
problem.internalJ = internalJ;
problem.depotEntryL = zeros(R, 2);
problem.depotEntryH = zeros(R, 2);
problem.depotEntryJ = zeros(R, 2);
problem.exitDepotL = zeros(R, 2);
problem.exitDepotH = zeros(R, 2);
problem.exitDepotJ = zeros(R, 2);
problem.transitionL = inf(R, 2, R, 2);
problem.transitionH = inf(R, 2, R, 2);
problem.transitionJ = inf(R, 2, R, 2);
for road = 1:R
    for direction = 1:2
        entryPoint = reshape(entry(road, direction, :), 1, 3);
        exitRoad = reshape(exitPoint(road, direction, :), 1, 3);
        [problem.depotEntryJ(road,direction), ...
            problem.depotEntryL(road, direction), ...
            problem.depotEntryH(road, direction)] = ...
            legMetrics(launchPoint,entryPoint,param);
        [problem.exitDepotJ(road,direction), ...
            problem.exitDepotL(road, direction), ...
            problem.exitDepotH(road, direction)] = ...
            legMetrics(exitRoad,recoveryPoint,param);
        for nextRoad = 1:R
            if nextRoad == road
                continue;
            end
            for nextDirection = 1:2
                nextEntry = reshape(entry( ...
                    nextRoad, nextDirection, :), 1, 3);
                [problem.transitionJ(road,direction,nextRoad,nextDirection), ...
                    problem.transitionL(road, direction, ...
                    nextRoad, nextDirection), ...
                    problem.transitionH(road, direction, ...
                    nextRoad, nextDirection)] = ...
                    legMetrics(exitRoad,nextEntry,param);
            end
        end
    end
end
end

function [energy_kJ,length3D,climb] = legMetrics(firstPoint,secondPoint,param)
[energy_kJ,length3D,climb] = ...
    uav_edge_energy(firstPoint,secondPoint,param.energy);
end

function solution = evaluateSolution(order, directions, problem)
R = problem.R;
firstRoad = order(1);
firstDirection = directions(1);
length3D = problem.depotEntryL(firstRoad, firstDirection);
climb = problem.depotEntryH(firstRoad, firstDirection);
energy = problem.depotEntryJ(firstRoad,firstDirection);
time_s = problem.depotEntryTime(firstRoad,firstDirection);
transferEnergy = energy;
transferLength = length3D;
detourCount = double(problem.depotEntryDetour(firstRoad,firstDirection));
unsafeDirectCount = double(problem.depotEntryDirectMinimumClearance( ...
    firstRoad,firstDirection)<problem.requiredClearance-1e-6);
minimumClearance = problem.depotEntryMinimumClearance(firstRoad,firstDirection);
for position = 1:R
    road = order(position);
    direction = directions(position);
    length3D = length3D+problem.internalL(road, direction);
    climb = climb+problem.internalH(road, direction);
    energy = energy+problem.internalJ(road,direction);
    time_s = time_s+problem.internalTime(road,direction);
    if position < R
        nextRoad = order(position+1);
        nextDirection = directions(position+1);
        length3D = length3D+problem.transitionL( ...
            road, direction, nextRoad, nextDirection);
        climb = climb+problem.transitionH( ...
            road, direction, nextRoad, nextDirection);
        energy = energy+problem.transitionJ( ...
            road,direction,nextRoad,nextDirection);
        time_s = time_s+problem.transitionTime( ...
            road,direction,nextRoad,nextDirection);
        transferEnergy = transferEnergy+problem.transitionJ( ...
            road,direction,nextRoad,nextDirection);
        transferLength = transferLength+problem.transitionL( ...
            road,direction,nextRoad,nextDirection);
        detourCount = detourCount+double(problem.transitionDetour( ...
            road,direction,nextRoad,nextDirection));
        unsafeDirectCount = unsafeDirectCount+double( ...
            problem.transitionDirectMinimumClearance( ...
            road,direction,nextRoad,nextDirection)< ...
            problem.requiredClearance-1e-6);
        minimumClearance = min(minimumClearance, ...
            problem.transitionMinimumClearance( ...
            road,direction,nextRoad,nextDirection));
    else
        length3D = length3D+problem.exitDepotL(road, direction);
        climb = climb+problem.exitDepotH(road, direction);
        energy = energy+problem.exitDepotJ(road,direction);
        time_s = time_s+problem.exitDepotTime(road,direction);
        transferEnergy = transferEnergy+problem.exitDepotJ(road,direction);
        transferLength = transferLength+problem.exitDepotL(road,direction);
        detourCount = detourCount+double(problem.exitDepotDetour(road,direction));
        unsafeDirectCount = unsafeDirectCount+double( ...
            problem.exitDepotDirectMinimumClearance(road,direction)< ...
            problem.requiredClearance-1e-6);
        minimumClearance = min(minimumClearance, ...
            problem.exitDepotMinimumClearance(road,direction));
    end
end
solution.order = double(order(:)');
solution.directions = double(directions(:)');
solution.L3D = length3D;
solution.Hplus = climb;
solution.J = energy;
solution.Time_s = time_s;
solution.TransferEnergy_kJ = transferEnergy;
solution.TransferLength_m = transferLength;
solution.DetourCount = detourCount;
solution.UnsafeDirectCount = unsafeDirectCount;
solution.MinimumTransferClearance_m = minimumClearance;
end

function cost = fixedOrderCost(order, problem)
firstRoad = order(1);
previous = problem.depotEntryJ(firstRoad, :)+ ...
    problem.internalJ(firstRoad, :);
for position = 2:problem.R
    road = order(position);
    precedingRoad = order(position-1);
    current = inf(1, 2);
    for direction = 1:2
        transition = [problem.transitionJ( ...
            precedingRoad, 1, road, direction), ...
            problem.transitionJ(precedingRoad, 2, road, direction)];
        current(direction) = min(previous+transition)+ ...
            problem.internalJ(road, direction);
    end
    previous = current;
end
cost = min(previous+problem.exitDepotJ(order(end), :));
end

function solution = solveFixedOrderDP(order, problem)
R = problem.R;
dp = inf(R, 2);
parent = zeros(R, 2, 'uint8');
firstRoad = order(1);
dp(1, :) = problem.depotEntryJ(firstRoad, :)+ ...
    problem.internalJ(firstRoad, :);
for position = 2:R
    road = order(position);
    precedingRoad = order(position-1);
    for direction = 1:2
        transition = [problem.transitionJ( ...
            precedingRoad, 1, road, direction), ...
            problem.transitionJ(precedingRoad, 2, road, direction)];
        [best, predecessor] = min(dp(position-1, :)+transition);
        dp(position, direction) = best+problem.internalJ(road, direction);
        parent(position, direction) = uint8(predecessor);
    end
end
[~, lastDirection] = min(dp(R, :)+problem.exitDepotJ(order(end), :));
directions = zeros(1, R);
directions(R) = lastDirection;
for position = R:-1:2
    directions(position-1) = parent(position, directions(position));
end
solution = evaluateSolution(order, directions, problem);
assert(costsEqual(solution.J, fixedOrderCost(order, problem)));
end

function order = nearestNeighbourOrder(problem)
remaining = true(1, problem.R);
order = zeros(1, problem.R);
currentRoad = 0;
currentDirection = 0;
for position = 1:problem.R
    bestCost = inf;
    bestRoad = 0;
    bestDirection = 0;
    for road = find(remaining)
        for direction = 1:2
            if position == 1
                cost = problem.depotEntryJ(road, direction);
            else
                cost = problem.transitionJ( ...
                    currentRoad, currentDirection, road, direction);
            end
            if cost < bestCost
                bestCost = cost;
                bestRoad = road;
                bestDirection = direction;
            end
        end
    end
    order(position) = bestRoad;
    remaining(bestRoad) = false;
    currentRoad = bestRoad;
    currentDirection = bestDirection;
end
end

function solution = solveGaDp(problem, ga, seed, heuristicSolution)
stream = RandStream('mt19937ar', 'Seed', seed);
R = problem.R;
populationSize = ga.populationSize;
population = zeros(populationSize, R);
population(1, :) = 1:R;
population(2, :) = heuristicSolution.order;
for individual = 3:populationSize
    population(individual, :) = randperm(stream, R);
end
fitness = evaluatePopulation(population, problem);
[bestCost, bestIndex] = min(fitness);
bestOrder = population(bestIndex, :);
stall = 0;
for generation = 1:ga.maximumGenerations
    [~, ranking] = sort(fitness, 'ascend');
    nextPopulation = zeros(size(population));
    elite = ranking(1:ga.eliteCount);
    nextPopulation(1:ga.eliteCount, :) = population(elite, :);
    writeIndex = ga.eliteCount+1;
    while writeIndex <= populationSize
        firstParent = tournamentSelect(fitness, ga.tournamentSize, stream);
        secondParent = tournamentSelect(fitness, ga.tournamentSize, stream);
        child1 = population(firstParent, :);
        child2 = population(secondParent, :);
        if rand(stream) < ga.crossoverRate
            [child1, child2] = orderCrossover(child1, child2, stream);
        end
        child1 = mutatePermutation( ...
            child1, ga.permutationMutationRate, stream);
        child2 = mutatePermutation( ...
            child2, ga.permutationMutationRate, stream);
        nextPopulation(writeIndex, :) = child1;
        if writeIndex+1 <= populationSize
            nextPopulation(writeIndex+1, :) = child2;
        end
        writeIndex = writeIndex+2;
    end
    population = nextPopulation;
    fitness = evaluatePopulation(population, problem);
    [generationBest, generationIndex] = min(fitness);
    if generationBest < bestCost- ...
            ga.improvementTolerance*max(1, abs(bestCost))
        bestCost = generationBest;
        bestOrder = population(generationIndex, :);
        stall = 0;
    else
        stall = stall+1;
    end
    if generation >= ga.minimumGenerations && stall >= ga.stallGenerations
        break;
    end
end
solution = solveFixedOrderDP(bestOrder, problem);
assert(costsEqual(solution.J, bestCost));
end

function fitness = evaluatePopulation(population, problem)
fitness = inf(size(population, 1), 1);
for individual = 1:size(population, 1)
    fitness(individual) = fixedOrderCost(population(individual, :), problem);
end
end

function index = tournamentSelect(fitness, tournamentSize, stream)
competitors = randi(stream, numel(fitness), 1, tournamentSize);
[~, localIndex] = min(fitness(competitors));
index = competitors(localIndex);
end

function [child1, child2] = orderCrossover(parent1, parent2, stream)
cuts = sort(randperm(stream, numel(parent1), 2));
child1 = oxChild(parent1, parent2, cuts(1), cuts(2));
child2 = oxChild(parent2, parent1, cuts(1), cuts(2));
end

function child = oxChild(primary, secondary, first, last)
R = numel(primary);
child = zeros(1, R);
child(first:last) = primary(first:last);
slots = [last+1:R, 1:first-1];
scan = [secondary(last+1:R), secondary(1:last)];
remaining = scan(~ismember(scan, child(first:last)));
child(slots) = remaining;
end

function chromosome = mutatePermutation(chromosome, rate, stream)
if rand(stream) >= rate
    return;
end
positions = sort(randperm(stream, numel(chromosome), 2));
if rand(stream) < 0.5
    chromosome(positions(1):positions(2)) = ...
        chromosome(positions(2):-1:positions(1));
else
    chromosome(positions) = chromosome(fliplr(positions));
end
end

function validateSolution(solution, problem)
assert(numel(solution.order) == problem.R && ...
    numel(unique(solution.order)) == problem.R);
assert(all(sort(solution.order) == 1:problem.R));
assert(all(ismember(solution.directions, [1, 2])));
assert(all(isfinite([solution.J,solution.L3D,solution.Hplus, ...
    solution.Time_s,solution.TransferEnergy_kJ, ...
    solution.MinimumTransferClearance_m])));
assert(solution.L3D > 0 && solution.Hplus >= 0);
assert(solution.MinimumTransferClearance_m >= ...
    problem.requiredClearance-1e-5);
recalculated = evaluateSolution( ...
    solution.order, solution.directions, problem);
assert(costsEqual(solution.J, recalculated.J));
end

function yes = costsEqual(first, second)
yes = abs(first-second) <= 1e-8*max(1, max(abs([first, second])));
end

%% Measured-road reconstruction and satellite crop

function path = findInput(dataDir, candidates)
for index = 1:numel(candidates)
    matches = dir(fullfile(dataDir, '**', candidates{index}));
    if ~isempty(matches)
        path = fullfile(matches(1).folder, matches(1).name);
        return;
    end
end
error('run_exp2b:MissingInput', ...
    'Missing input in %s. Expected one of: %s', ...
    dataDir, strjoin(candidates, ', '));
end

function centreline = reconstructCentreline(record, majorShapes, ...
    minorShapes, param)
if strcmp(record.sourceLayer, 'Major')
    shape = majorShapes(record.featureId);
else
    shape = minorShapes(record.featureId);
end
[x, y] = longestFinitePart(shape.X, shape.Y);
duplicate = [false; hypot(diff(x), diff(y)) < 1e-6];
x(duplicate) = [];
y(duplicate) = [];
sourceProgress = [0; cumsum(hypot(diff(x), diff(y)))];
progress = (0:param.centerlineSpacing:sourceProgress(end))';
if progress(end) < sourceProgress(end)
    progress(end+1, 1) = sourceProgress(end);
end
xq = interp1(sourceProgress, x, progress, 'pchip');
yq = interp1(sourceProgress, y, progress, 'pchip');
localProgress = (0:param.centerlineSpacing:param.routeLength)';
query = record.startOffset+localProgress;
assert(query(end) <= progress(end)+1e-8);
centreline.s = localProgress;
centreline.x = interp1(progress, xq, query, 'pchip');
centreline.y = interp1(progress, yq, query, 'pchip');
assert(all(isfinite([centreline.x; centreline.y])));
end

function [x, y] = longestFinitePart(X, Y)
finitePoints = isfinite(X) & isfinite(Y);
edges = diff([false, finitePoints, false]);
starts = find(edges == 1);
stops = find(edges == -1)-1;
assert(~isempty(starts));
[~, index] = max(stops-starts+1);
x = X(starts(index):stops(index));
y = Y(starts(index):stops(index));
x = x(:);
y = y(:);
end

function [xLimits, yLimits] = commonMapLimits(indices, ...
    centrelineByLibrary,launchPoint,recoveryPoint,reference,marginFraction)
points = [launchPoint(1:2);recoveryPoint(1:2)];
for libraryIndex = indices
    centreline = centrelineByLibrary{libraryIndex};
    points = [points; centreline.x, centreline.y]; %#ok<AGROW>
end
minimum = min(points, [], 1);
maximum = max(points, [], 1);
span = maximum-minimum;
span = max(span, max(span));
margin = max(40, marginFraction*max(span));
xLimits = [minimum(1)-margin, maximum(1)+margin];
yLimits = [minimum(2)-margin, maximum(2)+margin];
xLimits = fitLimitsToReference(xLimits, reference.XWorldLimits);
yLimits = fitLimitsToReference(yLimits, reference.YWorldLimits);
end

function limits = fitLimitsToReference(limits, available)
if diff(limits) > diff(available)
    limits = available;
    return;
end
if limits(1) < available(1)
    limits = limits+(available(1)-limits(1));
end
if limits(2) > available(2)
    limits = limits-(limits(2)-available(2));
end
end

function [crop, cropReference] = cropMap(image, reference, xLimits, yLimits)
dx = reference.CellExtentInWorldX;
dy = reference.CellExtentInWorldY;
column1 = max(1, floor((xLimits(1)-reference.XWorldLimits(1))/dx)+1);
column2 = min(size(image, 2), ...
    ceil((xLimits(2)-reference.XWorldLimits(1))/dx));
row1 = max(1, floor((reference.YWorldLimits(2)-yLimits(2))/dy)+1);
row2 = min(size(image, 1), ...
    ceil((reference.YWorldLimits(2)-yLimits(1))/dy));
crop = image(row1:row2, column1:column2, :);
xWorld = [reference.XWorldLimits(1)+(column1-1)*dx, ...
    reference.XWorldLimits(1)+column2*dx];
yWorld = [reference.YWorldLimits(2)-row2*dy, ...
    reference.YWorldLimits(2)-(row1-1)*dy];
cropReference = maprefcells(xWorld, yWorld, ...
    [size(crop, 1), size(crop, 2)], 'ColumnsStartFrom', 'north');
end

%% Drawing and export validation

function [targetHandle, transferHandle, detourHandle,launchHandle, ...
    recoveryHandle,directionHandle] = ...
    drawTrajectoryPanel(ax, task, routeLibrary, ...
    centrelineByLibrary,launchPoint,recoveryPoint,cfg)
solution = task.solution;
globalOrder = task.libraryIndices(solution.order);

% Draw the actual terrain-aware transfer polylines retained by Section 3.3.
problem = task.problem;
for position = 1:task.R
    road = solution.order(position);
    direction = solution.directions(position);
    if position == 1
        path = problem.depotEntryPath{road,direction};
        detour = problem.depotEntryDetour(road,direction);
    else
        precedingRoad = solution.order(position-1);
        precedingDirection = solution.directions(position-1);
        path = problem.transitionPath{precedingRoad,precedingDirection, ...
            road,direction};
        detour = problem.transitionDetour(precedingRoad,precedingDirection, ...
            road,direction);
    end
    if detour
        visibleDetour = drawDetourEvidence(ax,path,cfg);
    else
        visibleDetour = false;
    end
    drawTransferPath(ax,path,cfg,visibleDetour);
end
lastRoad = solution.order(end);
lastDirection = solution.directions(end);
lastPath = problem.exitDepotPath{lastRoad,lastDirection};
if problem.exitDepotDetour(lastRoad,lastDirection)
    visibleDetour = drawDetourEvidence(ax,lastPath,cfg);
else
    visibleDetour = false;
end
drawTransferPath(ax,lastPath,cfg,visibleDetour);

% Draw all target road centrelines in the same blue style as Figure 3.
for localRoad = 1:task.R
    libraryIndex = task.libraryIndices(localRoad);
    centreline = centrelineByLibrary{libraryIndex};
    plot(ax, centreline.x, centreline.y, '-', 'Color', 'w', ...
        'LineWidth', cfg.targetHaloWidth, 'HandleVisibility', 'off');
    plot(ax, centreline.x, centreline.y, '-', ...
        'Color', cfg.selectedColor, 'LineWidth', cfg.targetLineWidth, ...
        'HandleVisibility', 'off');
end

% One arrow per road shows its optimized direction.
for position = 1:task.R
    libraryIndex = globalOrder(position);
    centreline = centrelineByLibrary{libraryIndex};
    drawRoadDirectionArrow(ax, centreline, ...
        solution.directions(position), cfg);
end
% Visit labels retain the blue-box style of the road-network figure. They
    % are centred on the road; arrows are placed near the oriented entry.
drawVisitOrderLabels(ax, globalOrder, solution.directions, ...
    centrelineByLibrary,[launchPoint;recoveryPoint],cfg);

plot(ax,launchPoint(1),launchPoint(2),'p','MarkerSize',8.0, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.8, 'HandleVisibility', 'off');
plot(ax,launchPoint(1),launchPoint(2),'p','MarkerSize',6.2, ...
    'MarkerFaceColor',cfg.launchColor,'MarkerEdgeColor',[0.02,0.20,0.05], ...
    'LineWidth', 0.6, 'HandleVisibility', 'off');
plot(ax,recoveryPoint(1),recoveryPoint(2),'s','MarkerSize',7.2, ...
    'MarkerFaceColor','w','MarkerEdgeColor','w','LineWidth',1.8, ...
    'HandleVisibility','off');
plot(ax,recoveryPoint(1),recoveryPoint(2),'s','MarkerSize',5.7, ...
    'MarkerFaceColor',cfg.recoveryColor,'MarkerEdgeColor',[0.25,0.02,0.01], ...
    'LineWidth',0.6,'HandleVisibility','off');

% Dummy handles provide one shared legend without adding visible geometry.
targetHandle = plot(ax, nan, nan, '-', 'Color', cfg.selectedColor, ...
    'LineWidth', cfg.targetLineWidth);
    transferHandle = plot(ax, nan, nan, '--', 'Color', cfg.transferColor, ...
        'LineWidth', cfg.transferLineWidth);
    detourHandle = plot(ax,nan,nan,'-','Color',cfg.detourColor, ...
        'LineWidth',1.45);
launchHandle = plot(ax,nan,nan,'p','MarkerSize',6.2, ...
    'MarkerFaceColor',cfg.launchColor,'MarkerEdgeColor',[0.02,0.20,0.05]);
recoveryHandle = plot(ax,nan,nan,'s','MarkerSize',5.7, ...
    'MarkerFaceColor',cfg.recoveryColor,'MarkerEdgeColor',[0.25,0.02,0.01]);
directionHandle = plot(ax, nan, nan, '-', 'Color', cfg.directionColor, ...
    'LineWidth', 1.45, 'Marker', '>', 'MarkerSize', 6.5, ...
    'MarkerFaceColor', cfg.directionColor, 'MarkerEdgeColor', [0.25, 0.20, 0.00]);
end

function visibleDetour = drawDetourEvidence(ax,path,cfg)
plot(ax,path([1,end],1),path([1,end],2),':', ...
    'Color',cfg.unsafeChordColor,'LineWidth',0.8,'HandleVisibility','off');
a = path(1,1:2);
v = path(end,1:2)-a;
if norm(v) > 1e-9
    offset = abs((path(:,1)-a(1))*v(2)-(path(:,2)-a(2))*v(1))/norm(v);
else
    offset = vecnorm(path(:,1:2)-a,2,2);
end
visibleDetour = max(offset) >= cfg.visibleDetourThreshold_m;
end

function drawTransferPath(ax,path,cfg,visibleDetour)
plot(ax,path(:,1),path(:,2),'--','Color','w', ...
    'LineWidth',cfg.transferHaloWidth,'HandleVisibility','off');
if visibleDetour
    plot(ax,path(:,1),path(:,2),'-','Color',cfg.detourColor, ...
        'LineWidth',1.45,'HandleVisibility','off');
else
    plot(ax,path(:,1),path(:,2),'--','Color',cfg.transferColor, ...
        'LineWidth',cfg.transferLineWidth,'HandleVisibility','off');
end
end

function [entryPoint, exitPoint] = orientedEndpoints(route, direction)
if direction == 1
    entryPoint = route.waypoints(1, :);
    exitPoint = route.waypoints(end, :);
else
    entryPoint = route.waypoints(end, :);
    exitPoint = route.waypoints(1, :);
end
end

function drawRoadDirectionArrow(ax, centreline, direction, cfg)
if direction == 1
    xRoad = centreline.x(:);
    yRoad = centreline.y(:);
else
    xRoad = flipud(centreline.x(:));
    yRoad = flipud(centreline.y(:));
end
count = numel(xRoad);
halfSpan = max(5, round(0.055*count));
centre = max(halfSpan+1,round(0.18*count));
first = max(1, centre-halfSpan);
last = min(count, centre+halfSpan);
x = xRoad(first);
y = yRoad(first);
dx = xRoad(last)-x;
dy = yRoad(last)-y;
arrowLength = hypot(dx, dy);
if arrowLength <= eps
    return;
end
forward = [dx, dy]/arrowLength;
normal = [-forward(2), forward(1)];
tip = [x+dx, y+dy];
headLength = 0.42*arrowLength;
headHalfWidth = 0.25*arrowLength;
base = tip-headLength*forward;
whiteBase = tip-1.18*headLength*forward;
whiteHalfWidth = 1.30*headHalfWidth;

plot(ax, [x, base(1)], [y, base(2)], '-', 'Color', 'w', ...
    'LineWidth', 3.5, 'HandleVisibility', 'off');
patch(ax, [tip(1), whiteBase(1)+whiteHalfWidth*normal(1), ...
    whiteBase(1)-whiteHalfWidth*normal(1)], ...
    [tip(2), whiteBase(2)+whiteHalfWidth*normal(2), ...
    whiteBase(2)-whiteHalfWidth*normal(2)], 'w', ...
    'EdgeColor', 'w', 'HandleVisibility', 'off');
plot(ax, [x, base(1)], [y, base(2)], '-', ...
    'Color', cfg.directionColor, 'LineWidth', 1.65, ...
    'HandleVisibility', 'off');
patch(ax, [tip(1), base(1)+headHalfWidth*normal(1), ...
    base(1)-headHalfWidth*normal(1)], ...
    [tip(2), base(2)+headHalfWidth*normal(2), ...
    base(2)-headHalfWidth*normal(2)], cfg.directionColor, ...
    'EdgeColor', [0.25, 0.20, 0.00], 'LineWidth', 0.60, ...
    'HandleVisibility', 'off');
end

function drawVisitOrderLabels(ax, globalOrder, directions, ...
    centrelineByLibrary,boundaryPoints,cfg)
xLimits = xlim(ax);
yLimits = ylim(ax);
placed = zeros(numel(globalOrder), 3);
for visitPosition = 1:numel(globalOrder)
    centreline = centrelineByLibrary{globalOrder(visitPosition)};
    if directions(visitPosition) == 1
        x = centreline.x(:);
        y = centreline.y(:);
    else
        x = flipud(centreline.x(:));
        y = flipud(centreline.y(:));
    end
    candidates = visitLabelCandidates(x, y);
    halfWidth = 0.023+0.012*(visitPosition >= 10);
    [labelPoint, normalizedPoint] = chooseVisitLabelCandidate( ...
        candidates, halfWidth, placed(1:visitPosition-1, :), ...
        boundaryPoints,xLimits,yLimits,cfg.showMetrics);
    placed(visitPosition, :) = [normalizedPoint, halfWidth];
    text(ax, labelPoint(1), labelPoint(2), sprintf('%d', visitPosition), ...
        'Color', 'w', 'FontName', 'Times New Roman', 'FontWeight', 'bold', ...
        'FontSize', cfg.routeLabelFontSize, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'BackgroundColor', cfg.selectedColor, ...
        'Margin', 1.0, 'Clipping', 'on', 'HandleVisibility', 'off');
end
end

function candidates = visitLabelCandidates(x, y)
count = numel(x);
% Put the visit sequence at the road midpoint.  Nearby candidates are used
% only when the midpoint label would overlap another annotation.
fractions = [0.50, 0.44, 0.56, 0.38, 0.62, 0.32, 0.68];
candidates = zeros(2*numel(fractions), 2);
candidateIndex = 1;
for fraction = fractions
    centreIndex = max(2, min(count-1, round(1+fraction*(count-1))));
    before = max(1, centreIndex-2);
    after = min(count, centreIndex+2);
    tangent = [x(after)-x(before), y(after)-y(before)];
    tangent = tangent/max(norm(tangent), eps);
    normal = [-tangent(2), tangent(1)];
    point = [x(centreIndex), y(centreIndex)];
    candidates(candidateIndex, :) = point+16*normal;
    candidates(candidateIndex+1, :) = point-16*normal;
    candidateIndex = candidateIndex+2;
end
end

function [point, normalizedPoint] = chooseVisitLabelCandidate( ...
    candidates,halfWidth,placed,boundaryPoints,xLimits,yLimits,showMetrics)
normalized = [(candidates(:, 1)-xLimits(1))/diff(xLimits), ...
    (candidates(:, 2)-yLimits(1))/diff(yLimits)];
boundaryNormalized = [(boundaryPoints(:,1)-xLimits(1))/diff(xLimits), ...
    (boundaryPoints(:,2)-yLimits(1))/diff(yLimits)];
protected = normalized(:, 1) >= 0.04 & normalized(:, 1) <= 0.96 & ...
    normalized(:, 2) >= 0.04 & normalized(:, 2) <= 0.96;
if showMetrics
    protected = protected & ...
        ~(normalized(:, 1) >= 0.53 & normalized(:, 2) <= 0.23);
end
protected = protected & ...
    ~(normalized(:, 1) <= 0.30 & normalized(:, 2) >= 0.90);
for boundaryIndex=1:size(boundaryNormalized,1)
    protected=protected & vecnorm(normalized- ...
        boundaryNormalized(boundaryIndex,:),2,2)>=0.055;
end
valid = protected;
for index = 1:size(placed, 1)
    separated = abs(normalized(:, 1)-placed(index, 1)) >= ...
        halfWidth+placed(index, 3)+0.008 | ...
        abs(normalized(:, 2)-placed(index, 2)) >= 0.060;
    valid = valid & separated;
end
choice = find(valid, 1, 'first');
if isempty(choice)
    choice = find(protected, 1, 'first');
end
if isempty(choice)
    fallback = normalized(:, 1) >= 0.03 & normalized(:, 1) <= 0.97 & ...
        normalized(:, 2) >= 0.03 & normalized(:, 2) <= 0.97;
    fallback = fallback & ...
        ~(normalized(:, 1) <= 0.30 & normalized(:, 2) >= 0.90);
    for boundaryIndex=1:size(boundaryNormalized,1)
        fallback=fallback & vecnorm(normalized- ...
            boundaryNormalized(boundaryIndex,:),2,2)>=0.045;
    end
    choices = find(fallback);
    if ~isempty(choices)
        [~, localChoice] = max(normalized(choices, 2)- ...
            0.25*normalized(choices, 1));
        choice = choices(localChoice);
    end
end
assert(~isempty(choice), 'Could not place a visit-order label inside the map.');
point = candidates(choice, :);
normalizedPoint = normalized(choice, :);
end

function finalAudit(tasks,summary,routeLibrary,launchPoint,recoveryPoint, ...
    xLimits, yLimits, pngPath, pdfPath)
assert(isequal([tasks.R], [5, 10, 20]));
assert(all(arrayfun(@(task) numel(unique(task.libraryIndices)) == task.R, ...
    tasks)));
assert(all(isfinite([launchPoint,recoveryPoint])) && ...
    norm(launchPoint-recoveryPoint)>1 && all(isfinite(summary.J_E_kJ)) && ...
    all(isfinite(summary.L3D_m)) && all(isfinite(summary.Hplus_m)));
assert(all(summary.J_E_kJ > 0) && all(summary.L3D_m > 0) && ...
    all(summary.Hplus_m >= 0));
assert(diff(xLimits) > 0 && diff(yLimits) > 0);
assert(isfile(pngPath) && isfile(pdfPath));
pngInformation = dir(pngPath);
pdfInformation = dir(pdfPath);
assert(pngInformation.bytes > 1e5 && pdfInformation.bytes > 1e4);
for taskIndex = 1:numel(tasks)
    assert(all(tasks(taskIndex).libraryIndices >= 1 & ...
        tasks(taskIndex).libraryIndices <= numel(routeLibrary)));
end
fprintf('Final audit passed: representative sets, fixed open boundaries, finite metrics, ');
fprintf('measured-road provenance, and both figure exports are valid.\n');
end
