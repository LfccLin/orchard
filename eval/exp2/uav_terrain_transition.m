function result = uav_terrain_transition(startPoint,endPoint,startVelocity, ...
    endVelocity,terrain,param)
%UAV_TERRAIN_TRANSITION Minimum-energy terrain-safe position transfer.
% Coordinates are EPSG:param.crs metres and ellipsoidal/DEM heights in m.
% The road-end velocity arguments are retained for call compatibility but do
% not constrain the transfer.  The returned polyline is the exact layered
% graph path used for energy evaluation; its vertices are not smoothed away.

arguments
    startPoint (1,3) double {mustBeFinite}
    endPoint (1,3) double {mustBeFinite}
    startVelocity (1,3) double {mustBeFinite}
    endVelocity (1,3) double {mustBeFinite}
    terrain struct
    param struct
end

required = {'longitudinalStep','lateralStep','lateralHalfWidth', ...
    'verticalOffsets','edgeSampleStep','clearance','demErrorMargin'};
assert(all(isfield(param.terrain,required)), ...
    'Terrain-planning parameters are incomplete.');

directPath = [startPoint;endPoint];
[directSafe,directMinimumClearance] = pathIsSafe(directPath,terrain,param);
direct = evaluatePolyline(directPath,param,false);
direct.path = directPath;
direct.mode = "direct";
direct.minimumClearance_m = directMinimumClearance;

% Explicitly discard endpoint tangents.  Multirotor road-to-road transfer is
% represented as a position-to-position problem; transient reorientation at
% a graph vertex is outside this steady-segment energy model.
assert(all(isfinite([startVelocity,endVelocity])), ...
    'Compatibility velocity inputs must be finite even though unconstrained.');

horizontalDelta = endPoint(1:2)-startPoint(1:2);
horizontalDistance = norm(horizontalDelta);
if horizontalDistance < 1e-6
    fallbackPath = verticalFallback(startPoint,endPoint,terrain,param);
    graphCandidate = evaluatePolyline(fallbackPath,param,false);
    graphCandidate.path = fallbackPath;
    graphCandidate.mode = "raised-overflight";
    [~,graphCandidate.minimumClearance_m] = pathIsSafe( ...
        fallbackPath,terrain,param);
else
    graphCandidate = layeredShortestPath(startPoint,endPoint,terrain,param);
end

raisedPath = verticalFallback(startPoint,endPoint,terrain,param);
raisedCandidate = evaluatePolyline(raisedPath,param,false);
raisedCandidate.path = raisedPath;
raisedCandidate.mode = "raised-overflight";
[raisedSafe,raisedCandidate.minimumClearance_m] = pathIsSafe( ...
    raisedPath,terrain,param);
assert(raisedSafe,'The deterministic raised-overflight candidate is unsafe.');

if directSafe
    candidateSet = {graphCandidate,raisedCandidate,direct};
else
    candidateSet = {graphCandidate,raisedCandidate};
end
candidateEnergy = cellfun(@(candidate) candidate.energy_kJ,candidateSet);
[~,bestCandidate] = min(candidateEnergy);
result = candidateSet{bestCandidate};

result.directSafe = directSafe;
result.directMinimumClearance_m = directMinimumClearance;
result.detour = result.mode ~= "direct";
result.lateralDeviation_m = maximumPlanarDeviation(result.path);
result.turnCount = countPlanarTurns(result.path);
if ~result.detour
    result.detourMode = "none";
elseif result.lateralDeviation_m >= param.terrain.lateralStep-1e-6
    result.detourMode = "lateral-bypass";
else
    result.detourMode = "vertical-overflight";
end
assert(isfinite(result.energy_kJ) && result.energy_kJ >= 0);
assert(result.minimumClearance_m >= ...
    param.terrain.clearance+param.terrain.demErrorMargin-1e-5);
end

function result = layeredShortestPath(p0,p1,terrain,param)
delta = p1(1:2)-p0(1:2);
D = norm(delta);
forward = delta/D;
normal = [-forward(2),forward(1)];
layerCount = max(2,ceil(D/param.terrain.longitudinalStep));
xi = linspace(0,D,layerCount+1);

halfWidth = min(param.terrain.lateralHalfWidth, ...
    max(2*param.terrain.lateralStep,0.30*D));
lateralCount = floor(halfWidth/param.terrain.lateralStep);
eta = (-lateralCount:lateralCount)*param.terrain.lateralStep;
zeta = param.terrain.verticalOffsets(:)';

layers = cell(layerCount+1,1);
layers{1} = p0;
layers{end} = p1;
requiredClearance = param.terrain.clearance+param.terrain.demErrorMargin;
for g = 2:layerCount
    xy = p0(1:2)+xi(g)*forward+eta(:)*normal;
    ground = terrainHeight(xy,terrain,param.crs);
    valid = isfinite(ground);
    xy = xy(valid,:);
    ground = ground(valid);
    nodes = zeros(numel(ground)*numel(zeta),3);
    cursor = 0;
    for m = 1:numel(ground)
        for k = 1:numel(zeta)
            cursor = cursor+1;
            nodes(cursor,:) = [xy(m,:),ground(m)+requiredClearance+zeta(k)];
        end
    end
    layers{g} = nodes;
end

cost = {0};
parent = cell(layerCount+1,1);
for g = 2:layerCount+1
    current = layers{g};
    previous = layers{g-1};
    currentCost = inf(size(current,1),1);
    currentParent = zeros(size(current,1),1,'uint16');
    for v = 1:size(current,1)
        for u = 1:size(previous,1)
            if ~isfinite(cost{g-1}(u))
                continue;
            end
            lateralChange = abs(dot(current(v,1:2)-previous(u,1:2),normal));
            if lateralChange > 2.01*param.terrain.lateralStep
                continue;
            end
            edge = [previous(u,:);current(v,:)];
            [safe,~] = pathIsSafe(edge,terrain,param);
            if ~safe
                continue;
            end
            metric = evaluatePolyline(edge,param,false);
            trial = cost{g-1}(u)+metric.energy_kJ;
            if trial < currentCost(v)
                currentCost(v) = trial;
                currentParent(v) = uint16(u);
            end
        end
    end
    cost{g} = currentCost;
    parent{g} = currentParent;
end

[best,last] = min(cost{end});
if ~isfinite(best)
    result = evaluatePolyline(verticalFallback(p0,p1,terrain,param),param,false);
    result.path = verticalFallback(p0,p1,terrain,param);
    result.mode = "raised-fallback";
    [~,result.minimumClearance_m] = pathIsSafe(result.path,terrain,param);
    return;
end

indices = zeros(layerCount+1,1);
indices(end) = last;
for g = layerCount+1:-1:2
    indices(g-1) = parent{g}(indices(g));
end
path = zeros(layerCount+1,3);
for g = 1:layerCount+1
    path(g,:) = layers{g}(indices(g),:);
end
% Retain the selected node at every progress layer.  These vertices expose
% the actual graph exploration and make each terrain-driven heading change
% auditable in the trajectory figure.
path = removeNearDuplicates(path);
result = evaluatePolyline(path,param,false);
result.path = path;
result.mode = "layered-energy-DP";
[~,result.minimumClearance_m] = pathIsSafe(path,terrain,param);
end

function path = verticalFallback(p0,p1,terrain,param)
sampleCount = max(2,ceil(norm(p1(1:2)-p0(1:2))/ ...
    param.terrain.edgeSampleStep)+1);
t = linspace(0,1,sampleCount)';
xy = p0(1:2)+(p1(1:2)-p0(1:2)).*t;
ground = terrainHeight(xy,terrain,param.crs);
required = max(ground,[],'omitnan')+param.terrain.clearance+ ...
    param.terrain.demErrorMargin;
apex = max([p0(3),p1(3),required]);
path = [p0;p0(1:2),apex;p1(1:2),apex;p1];
path = removeCollinear(path);
end

function metric = evaluatePolyline(path,param,fixedPhotographySpeed)
energy = 0;
length3D = 0;
climb = 0;
time = 0;
for k = 1:size(path,1)-1
    if fixedPhotographySpeed
        [e,L,H,d] = uav_edge_energy(path(k,:),path(k+1,:), ...
            param.energy,'HorizontalSpeed',param.photographySpeed);
    else
        [e,L,H,d] = uav_edge_energy(path(k,:),path(k+1,:),param.energy);
    end
    energy = energy+e;
    length3D = length3D+L;
    climb = climb+H;
    time = time+d.time_s;
end
metric = struct('energy_kJ',energy,'length3D_m',length3D, ...
    'climb_m',climb,'time_s',time,'minimumClearance_m',NaN);
end

function [safe,minimumClearance] = pathIsSafe(path,terrain,param)
required = param.terrain.clearance+param.terrain.demErrorMargin;
minimumClearance = inf;
safe = true;
for k = 1:size(path,1)-1
    count = max(2,ceil(norm(path(k+1,1:2)-path(k,1:2))/ ...
        param.terrain.edgeSampleStep)+1);
    t = linspace(0,1,count)';
    xyz = path(k,:)+(path(k+1,:)-path(k,:)).*t;
    ground = terrainHeight(xyz(:,1:2),terrain,param.crs);
    clearance = xyz(:,3)-ground;
    if any(~isfinite(clearance))
        safe = false;
        minimumClearance = -inf;
        return;
    end
    minimumClearance = min(minimumClearance,min(clearance));
    if any(clearance < required-1e-6)
        safe = false;
    end
end
end

function height = terrainHeight(xy,terrain,crsCode)
if isfield(terrain,'interpolant')
    height = terrain.interpolant(xy(:,1),xy(:,2));
    return;
end
[latitude,longitude] = projinv(projcrs(crsCode),xy(:,1),xy(:,2));
[xi,yi] = geographicToIntrinsic(terrain.reference,latitude,longitude);
height = interp2(terrain.dem,xi,yi,'linear',NaN);
end

function path = removeCollinear(path)
path = removeNearDuplicates(path);
if size(path,1) <= 2
    return;
end
keep = true(size(path,1),1);
for k = 2:size(path,1)-1
    a = path(k,:)-path(k-1,:);
    b = path(k+1,:)-path(k,:);
    if norm(cross(a,b)) <= 1e-8*max(1,norm(a)*norm(b))
        keep(k) = false;
    end
end
path = path(keep,:);
end

function path = removeNearDuplicates(path)
if size(path,1) <= 1
    return;
end
path = path([true;vecnorm(diff(path,1,1),2,2)>1e-7],:);
end

function deviation = maximumPlanarDeviation(path)
xy = path(:,1:2);
origin = xy(1,:);
chord = xy(end,:)-origin;
if norm(chord) <= 1e-9
    deviation = max(vecnorm(xy-origin,2,2));
    return;
end
unit = chord/norm(chord);
projection = origin+((xy-origin)*unit').*unit;
deviation = max(vecnorm(xy-projection,2,2));
end

function count = countPlanarTurns(path)
if size(path,1) < 3
    count = 0;
    return;
end
segments = diff(path(:,1:2),1,1);
lengths = vecnorm(segments,2,2);
valid = lengths > 1e-7;
segments = segments(valid,:)./lengths(valid);
if size(segments,1) < 2
    count = 0;
    return;
end
cosine = sum(segments(1:end-1,:).*segments(2:end,:),2);
turnAngle = acos(max(-1,min(1,cosine)));
count = nnz(turnAngle >= deg2rad(5));
end
