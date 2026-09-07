function run_exp1b()
%RUN_EXP1B Experiment 1b: road geometry versus adaptive feasible spacing.
%
% This experiment does not modify or redraw the 3-D comparison figures from
% run_exp1. It reuses the road models, parameters, and proposed-adaptive
% results saved by run_exp1 in:
%   <project>/results/exp1/exp1_results.mat
%
% Outputs are written to the same directory:
%   <project>/results/exp1
%
% Main outputs:
%   fig_exp1b_local_geometry.png / .pdf
%   fig_exp1b_feasible_spacing.png / .pdf
%   exp1b_candidate_metrics.csv
%   exp1b_spearman.csv
%   exp1b_results.mat
%
% The analysis quantity is the candidate-wise farthest feasible forward
% spacing D*(s), already computed in run_exp1 with the road-aligned coverage
% bank and the same overlap / interval-coverage constraints:
%
%   D*(s) = max{d <= D0 : rho(s,s+d) >= rho0,
%                         C(s,s+d) >= Cmin }.
%
% Local geometry is evaluated in the same fixed forward D0 window:
%   K(s) = integral |d psi / ds| ds  -> discrete accumulated heading change
%   G(s) = integral |d beta / ds| ds -> discrete accumulated slope change
%
% The last D0 metres are excluded from D*, K, G correlation analysis so
% that road-end truncation is not mistaken for geometry-induced shortening.

close all;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
outputDir = fullfile(projectRoot, 'results', 'exp1');
exp1File = fullfile(outputDir, 'exp1_results.mat');

if ~isfile(exp1File)
    error(['Missing %s. Run run_exp1 first so that the road models, ', ...
        'parameters, and adaptive candidate-wise spacing are available.'], ...
        exp1File);
end

inputData = load(exp1File, 'roadModels', 'methods', 'param', 'results');
roadModels = inputData.roadModels;
methods = inputData.methods; %#ok<NASGU>
param = inputData.param;
results = inputData.results;

assert(numel(roadModels) == 3, 'Expected roads A, B, and C.');
assert(abs(param.centerlineSpacing-2) < 1e-10);
assert(abs(param.maximumSpacing-8.99) < 0.02);

%% Locate the proposed adaptive method saved by Exp. 1
methodIds = string({results.methodId});
proposedMethodIndex = find(methodIds == "c", 1);
if isempty(proposedMethodIndex)
    proposedMethodIndex = find(contains(lower(string({results.methodName})), ...
        "adaptive"), 1);
end
assert(~isempty(proposedMethodIndex), ...
    'Could not locate the proposed adaptive method in exp1_results.mat.');

proposedRoadResults = results(proposedMethodIndex).roads;

%% Candidate-wise geometry and adaptive-spacing diagnostics
analysis = repmat(struct( ...
    'roadId', '', ...
    'progress_m', [], ...
    'relativeElevation_m', [], ...
    'K_rad', [], ...
    'K_deg', [], ...
    'G_rad', [], ...
    'G_deg', [], ...
    'Dstar_m', [], ...
    'validMask', [], ...
    'adaptiveWaypointMask', [], ...
    'rho_D_K', NaN, ...
    'p_D_K', NaN, ...
    'rho_D_G', NaN, ...
    'p_D_G', NaN), 1, numel(roadModels));

candidateRows = cell(0, 9);
correlationRows = cell(0, 7);

for roadIndex = 1:numel(roadModels)
    model = roadModels(roadIndex);
    proposed = proposedRoadResults(roadIndex);

    if ~isfield(proposed, 'localMaximumFeasibleSpacing')
        error(['The saved Exp. 1 result does not contain ', ...
            'localMaximumFeasibleSpacing. Re-run the current run_exp1.m ', ...
            'once, then run run_exp1b again.']);
    end

    s = model.cumulativeDistance(:);
    zRelative = model.worldCoordinates(:, 3) - ...
        min(model.worldCoordinates(:, 3));

    % D*(s) comes from exactly the same road-aligned feasibility test used
    % by the proposed adaptive method in run_exp1.
    Dstar = proposed.localMaximumFeasibleSpacing(:);

    % Fixed forward D0-window geometry.
    K = localAccumulatedVariation(s, model.heading(:), param.maximumSpacing);
    G = localAccumulatedVariation(s, model.longitudinalSlope(:), ...
        param.maximumSpacing);

    % Exclude the final D0 window from mechanism analysis. A short D* there
    % can be caused only by the road ending, not by road geometry.
    validWindow = s <= s(end)-param.maximumSpacing + distanceTolerance(s);
    valid = validWindow & isfinite(Dstar) & Dstar > 0 & ...
        isfinite(K) & isfinite(G);

    if any(validWindow & (Dstar <= 0 | ~isfinite(Dstar)))
        warning('run_exp1b:NoFeasibleForwardCandidate', ...
            ['Road %s contains interior candidate positions with no ', ...
            'positive feasible forward spacing. They are excluded from ', ...
            'the correlation analysis.'], model.routeId);
    end

    adaptiveWaypointMask = false(size(s));
    adaptiveWaypointMask(proposed.selectedIndices) = true;

    [rhoDK, pDK] = spearmanWithP(Dstar(valid), K(valid));
    [rhoDG, pDG] = spearmanWithP(Dstar(valid), G(valid));

    analysis(roadIndex).roadId = model.routeId;
    analysis(roadIndex).progress_m = s;
    analysis(roadIndex).relativeElevation_m = zRelative;
    analysis(roadIndex).K_rad = K;
    analysis(roadIndex).K_deg = rad2deg(K);
    analysis(roadIndex).G_rad = G;
    analysis(roadIndex).G_deg = rad2deg(G);
    analysis(roadIndex).Dstar_m = Dstar;
    analysis(roadIndex).validMask = valid;
    analysis(roadIndex).adaptiveWaypointMask = adaptiveWaypointMask;
    analysis(roadIndex).rho_D_K = rhoDK;
    analysis(roadIndex).p_D_K = pDK;
    analysis(roadIndex).rho_D_G = rhoDG;
    analysis(roadIndex).p_D_G = pDG;

    for i = 1:numel(s)
        candidateRows(end+1, :) = { ...
            string(model.routeId), s(i), zRelative(i), ...
            rad2deg(K(i)), rad2deg(G(i)), Dstar(i), ...
            valid(i), adaptiveWaypointMask(i), ...
            remainingDistanceClass(s(i), s(end), param.maximumSpacing)}; ...
            %#ok<AGROW>
    end

    correlationRows(end+1, :) = { ...
        string(model.routeId), nnz(valid), rhoDK, pDK, rhoDG, pDG, ...
        mean(Dstar(valid), 'omitnan')}; %#ok<AGROW>
end

%% Pooled correlation (reported separately; road-wise values remain primary)
allD = [];
allK = [];
allG = [];
for roadIndex = 1:numel(analysis)
    valid = analysis(roadIndex).validMask;
    allD = [allD; analysis(roadIndex).Dstar_m(valid)]; %#ok<AGROW>
    allK = [allK; analysis(roadIndex).K_rad(valid)]; %#ok<AGROW>
    allG = [allG; analysis(roadIndex).G_rad(valid)]; %#ok<AGROW>
end
[pooledRhoDK, pooledPDK] = spearmanWithP(allD, allK);
[pooledRhoDG, pooledPDG] = spearmanWithP(allD, allG);
correlationRows(end+1, :) = { ...
    "Pooled", numel(allD), pooledRhoDK, pooledPDK, pooledRhoDG, pooledPDG, ...
    mean(allD, 'omitnan')};

candidateMetrics = cell2table(candidateRows, 'VariableNames', { ...
    'Road', 'Progress_m', 'RelativeElevation_m', ...
    'LocalTurning_K_deg', 'LocalSlopeVariation_G_deg', ...
    'Dstar_m', 'ValidForCorrelation', 'IsAdaptiveWaypoint', ...
    'RoadEndWindowClass'});

spearmanSummary = cell2table(correlationRows, 'VariableNames', { ...
    'Road', 'ValidCandidateCount', ...
    'Spearman_Dstar_K', 'p_Dstar_K', ...
    'Spearman_Dstar_G', 'p_Dstar_G', ...
    'MeanDstar_m'});

%% Publication figures: local geometry and feasible spacing saved separately
drawLocalGeometryFigure(analysis, outputDir);
drawFeasibleSpacingFigure(analysis, param, outputDir);

%% Save numerical outputs
writetable(candidateMetrics, ...
    fullfile(outputDir, 'exp1b_candidate_metrics.csv'));
writetable(spearmanSummary, ...
    fullfile(outputDir, 'exp1b_spearman.csv'));
save(fullfile(outputDir, 'exp1b_results.mat'), ...
    'analysis', 'candidateMetrics', 'spearmanSummary', 'param');

fprintf('Experiment 1b completed. Results written to:\n  %s\n', outputDir);
disp(spearmanSummary);
end


function accumulated = localAccumulatedVariation(s, signal, windowLength)
% Accumulated absolute signal change over [s_i, s_i + windowLength].
%
% For the sampled road, this is the discrete counterpart of:
%   integral |d signal / d xi| d xi
%
% The endpoint is linearly interpolated when s_i + windowLength falls
% between two samples, so K(s) and G(s) use the same physical D0 window
% rather than an index-count approximation.

s = s(:);
signal = signal(:);
n = numel(s);
accumulated = nan(n, 1);

assert(numel(signal) == n);
assert(all(diff(s) > 0));

for i = 1:n
    windowEnd = s(i) + windowLength;
    if windowEnd > s(end) + distanceTolerance(s)
        continue;
    end

    % Samples at or before the physical window end.
    lastFull = find(s <= windowEnd, 1, 'last');
    localS = s(i:lastFull);
    localSignal = signal(i:lastFull);

    if localS(end) < windowEnd - distanceTolerance(s)
        nextIndex = lastFull + 1;
        assert(nextIndex <= n);
        endpointValue = interp1( ...
            s([lastFull, nextIndex]), signal([lastFull, nextIndex]), ...
            windowEnd, 'linear');
        localSignal(end+1, 1) = endpointValue;
    end

    accumulated(i) = sum(abs(diff(localSignal)));
end
end


function [rho, p] = spearmanWithP(x, y)
x = x(:);
y = y(:);
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

if numel(x) < 3 || numel(unique(x)) < 2 || numel(unique(y)) < 2
    rho = NaN;
    p = NaN;
    return;
end

[rho, p] = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');
end


function label = remainingDistanceClass(s, roadEnd, D0)
if s <= roadEnd-D0 + distanceTolerance([0; roadEnd])
    label = "FullD0Window";
else
    label = "RoadEndExcluded";
end
end


function drawLocalGeometryFigure(analysis, outputDir)
% Figure 1: K(s) and G(s) for Roads A-C.

style.figureSize = [40, 80, 2500, 780];
style.fontSize = 20;
style.labelFontSize = 22;
style.panelFontSize = 23;
style.legendFontSize = 18;
style.lineWidth = 2.1;

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', style.figureSize);
set(fig, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman');

layout = tiledlayout(fig, 1, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

panelLetters = {'(a)', '(b)', '(c)'};

for roadIndex = 1:numel(analysis)
    a = analysis(roadIndex);
    s = a.progress_m;

    ax = nexttile(layout, roadIndex);
    hold(ax, 'on');

    yyaxis(ax, 'left');
    hK = plot(ax, s, a.K_deg, '-', ...
        'LineWidth', style.lineWidth);
    ylabel(ax, 'K(s) (deg)', ...
        'FontName', 'Times New Roman', ...
        'FontSize', style.labelFontSize);

    yyaxis(ax, 'right');
    hG = plot(ax, s, a.G_deg, '--', ...
        'LineWidth', style.lineWidth);
    ylabel(ax, 'G(s) (deg)', ...
        'FontName', 'Times New Roman', ...
        'FontSize', style.labelFontSize);

    xlabel(ax, 'Road progress s (m)', ...
        'FontName', 'Times New Roman', ...
        'FontSize', style.labelFontSize);

    xlim(ax, [s(1), s(end)]);
    grid(ax, 'on');
    box(ax, 'on');
    setCommonAxesStyle(ax, style.fontSize);

    % One panel label in the upper-left corner.
    text(ax, 0.018, 0.965, ...
        sprintf('%s Road %s', panelLetters{roadIndex}, a.roadId), ...
        'Units', 'normalized', ...
        'FontName', 'Times New Roman', ...
        'FontSize', style.panelFontSize, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left');

    % The legend meaning is identical for all roads, so show it once in a
    % dedicated row above the whole tiled figure to avoid covering data.
    if roadIndex == 2
        lgd = legend(ax, [hK, hG], ...
            {'Local turning K(s)', 'Slope variation G(s)'}, ...
            'Orientation', 'horizontal', ...
            'Box', 'on');
        lgd.Layout.Tile = 'north';
        set(lgd, 'FontName', 'Times New Roman', ...
            'FontSize', style.legendFontSize, ...
            'LineWidth', 1.0);
    end
end

fileStem = 'fig_exp1b_local_geometry';
exportgraphics(fig, fullfile(outputDir, [fileStem, '.png']), ...
    'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, fullfile(outputDir, [fileStem, '.pdf']), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');
close(fig);
end


function drawFeasibleSpacingFigure(analysis, param, outputDir)
% Figure 2: candidate-wise farthest feasible spacing D*(s) for Roads A-C.

style.figureSize = [40, 80, 2500, 780];
style.fontSize = 20;
style.labelFontSize = 22;
style.panelFontSize = 23;
style.legendFontSize = 18;
style.lineWidth = 2.1;
style.markerSize = 34;

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', style.figureSize);
set(fig, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman');

layout = tiledlayout(fig, 1, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

panelLetters = {'(a)', '(b)', '(c)'};

for roadIndex = 1:numel(analysis)
    a = analysis(roadIndex);
    s = a.progress_m;
    valid = a.validMask;
    waypointMask = a.adaptiveWaypointMask;

    ax = nexttile(layout, roadIndex);
    hold(ax, 'on');

    Dplot = a.Dstar_m;
    Dplot(~valid) = NaN;

    hD = stairs(ax, s, Dplot, '-', ...
        'LineWidth', style.lineWidth);

    % Reference spacing is shown as a dashed horizontal line. Its label is
    % written directly at the right side of the line rather than included
    % in the figure-wide legend.
    yline(ax, param.maximumSpacing, '--', ...
        'LineWidth', 1.8);

    hWp = scatter(ax, s(waypointMask & valid), ...
        Dplot(waypointMask & valid), ...
        style.markerSize, 'filled');

    xlim(ax, [s(1), s(end)]);
    ylim(ax, [0, max(10, ceil(param.maximumSpacing)+0.5)]);
    yticks(ax, [0, 2, 4, 6, 8]);
    grid(ax, 'on');
    box(ax, 'on');

    xlabel(ax, 'Road progress s (m)', ...
        'FontName', 'Times New Roman', ...
        'FontSize', style.labelFontSize);
    ylabel(ax, 'D^*(s) (m)', ...
        'FontName', 'Times New Roman', ...
        'FontSize', style.labelFontSize);

    setCommonAxesStyle(ax, style.fontSize);

    % Put the reference-spacing label at the right end of the dashed line.
    % A small inward x-offset keeps the text inside the axes boundary.
    xReferenceLabel = s(end) - 0.018*(s(end)-s(1));
    text(ax, xReferenceLabel, param.maximumSpacing, ...
        sprintf('Reference spacing D_0 = %.2f m', param.maximumSpacing), ...
        'Interpreter', 'tex', ...
        'FontName', 'Times New Roman', ...
        'FontSize', style.legendFontSize, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', ...
        'BackgroundColor', 'w', ...
        'Margin', 2, ...
        'Clipping', 'on');

    % One panel label in the upper-left corner.
    text(ax, 0.018, 0.965, ...
        sprintf('%s Road %s', panelLetters{roadIndex}, a.roadId), ...
        'Units', 'normalized', ...
        'FontName', 'Times New Roman', ...
        'FontSize', style.panelFontSize, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left');

    correlationText = sprintf( ...
        '\\rho_s(D^*,K)=%.2f,  \\rho_s(D^*,G)=%.2f', ...
        a.rho_D_K, a.rho_D_G);
    text(ax, 0.025, 0.075, correlationText, ...
        'Units', 'normalized', ...
        'Interpreter', 'tex', ...
        'FontName', 'Times New Roman', ...
        'FontSize', style.fontSize-1, ...
        'BackgroundColor', 'w', ...
        'EdgeColor', [0.25, 0.25, 0.25], ...
        'Margin', 4, ...
        'VerticalAlignment', 'bottom');

    % The legend meaning is identical for all roads, so show it once in a
    % dedicated row above the whole tiled figure. Reference spacing is not
    % included here because it is labelled directly beside the dashed line.
    if roadIndex == 2
        lgd = legend(ax, [hD, hWp], ...
            {'Maximum feasible spacing D^*(s)', ...
             'Adaptive waypoint'}, ...
            'Orientation', 'horizontal', ...
            'NumColumns', 2, ...
            'Box', 'on', ...
            'Interpreter', 'tex');
        lgd.Layout.Tile = 'north';
        set(lgd, 'FontName', 'Times New Roman', ...
            'FontSize', style.legendFontSize, ...
            'LineWidth', 1.0);
    end
end

fileStem = 'fig_exp1b_feasible_spacing';
exportgraphics(fig, fullfile(outputDir, [fileStem, '.png']), ...
    'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, fullfile(outputDir, [fileStem, '.pdf']), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');
close(fig);
end


function setCommonAxesStyle(ax, fontSize)
set(ax, 'FontName', 'Times New Roman', ...
    'FontSize', fontSize, ...
    'LineWidth', 1.0, ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'GridAlpha', 0.22);
end


function tolerance = distanceTolerance(s)
if isstruct(s)
    scale = max(1, s.cumulativeDistance(end));
else
    s = s(:);
    finiteS = s(isfinite(s));
    if isempty(finiteS)
        scale = 1;
    else
        scale = max(1, max(abs(finiteS)));
    end
end
tolerance = 128*eps(scale);
end
