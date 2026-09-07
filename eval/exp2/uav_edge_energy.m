function [energy_kJ, length3D_m, climb_m, details] = ...
    uav_edge_energy(firstPoint, secondPoint, energy, varargin)
%UAV_EDGE_ENERGY Directed steady-flight energy for one straight 3-D leg.
% The power model follows Gong et al., IEEE Transactions on Aerospace and
% Electronic Systems, 59(6), 7409-7422, 2023 (DOI: 10.1109/TAES.2023.3288846).
% Parameters must be supplied explicitly; this function invents no vehicle
% or terrain data. Output energy is kJ. By default, leg duration is chosen
% to minimize energy over the admissible speed interval. Pass
% 'HorizontalSpeed',v to impose a fixed horizontal imaging speed.

validateattributes(firstPoint, {'numeric'}, {'real','finite','vector','numel',3});
validateattributes(secondPoint, {'numeric'}, {'real','finite','vector','numel',3});
required = {'rotorCount','profileDrag','airDensity','rotorSolidity', ...
    'rotorDiscArea','thrustCoefficient','inducedCorrection','weight', ...
    'inducedVelocity','flatPlateHorizontal','flatPlateVertical', ...
    'horizontalSpeedMax','ascentSpeedMax','descentSpeedMax'};
assert(all(isfield(energy, required)), ...
    'The UAV energy parameter structure is incomplete.');

parser = inputParser;
parser.addParameter('HorizontalSpeed', NaN, ...
    @(x) isscalar(x) && (isnan(x) || (isfinite(x) && x > 0)));
parser.parse(varargin{:});
fixedHorizontalSpeed = parser.Results.HorizontalSpeed;
if ~isfield(energy, 'minimumSpeed')
    energy.minimumSpeed = 1;
end

delta = double(secondPoint(:).'-firstPoint(:).');
horizontal = hypot(delta(1), delta(2));
vertical = delta(3);
length3D_m = norm(delta);
climb_m = max(vertical, 0);
if length3D_m <= 1e-12
    energy_kJ = 0;
    details = struct('time_s',0,'horizontalSpeed_mps',0, ...
        'verticalSpeed_mps',0,'power_W',0);
    return;
end

minimumTime = max([horizontal/energy.horizontalSpeedMax, ...
    max(vertical,0)/energy.ascentSpeedMax, ...
    max(-vertical,0)/energy.descentSpeedMax]);

if isfinite(fixedHorizontalSpeed) && horizontal > 1e-9
    time_s = horizontal/fixedHorizontalSpeed;
    assert(time_s+1e-10 >= minimumTime, ...
        'Fixed imaging speed violates an ascent/descent speed limit.');
elseif horizontal > 1e-9
    slope = vertical/horizontal;
    horizontalUpper = energy.horizontalSpeedMax;
    if slope > 0
        horizontalUpper = min(horizontalUpper,energy.ascentSpeedMax/slope);
    elseif slope < 0
        horizontalUpper = min(horizontalUpper,energy.descentSpeedMax/(-slope));
    end
    horizontalLower = min(energy.minimumSpeed,horizontalUpper);
    horizontalSpeed = optimalHorizontalSpeed( ...
        slope,horizontalLower,horizontalUpper,energy);
    time_s = horizontal/horizontalSpeed;
else
    verticalUpper = energy.ascentSpeedMax;
    if vertical < 0
        verticalUpper = energy.descentSpeedMax;
    end
    verticalLower = min(energy.minimumSpeed,verticalUpper);
    verticalMagnitude = optimalVerticalSpeed( ...
        sign(vertical),verticalLower,verticalUpper,energy);
    time_s = abs(vertical)/verticalMagnitude;
end
horizontalSpeed = horizontal/time_s;
verticalSpeed = vertical/time_s;

power_W = powerAtDuration(time_s,horizontal,vertical,energy);
assert(isfinite(power_W) && power_W > 0 && time_s > 0);
energy_kJ = power_W*time_s/1000;
details = struct('time_s',time_s, ...
    'horizontalSpeed_mps',horizontalSpeed, ...
    'verticalSpeed_mps',verticalSpeed,'power_W',power_W);
end

function speed = optimalHorizontalSpeed(slope,lower,upper,energy)
persistent speedCache
if isempty(speedCache)
    speedCache = containers.Map('KeyType','char','ValueType','double');
end
key = sprintf('h|%.3f|%.3f|%.3f',slope,lower,upper);
if isKey(speedCache,key)
    speed = speedCache(key);
    return;
end
if upper <= lower*(1+1e-9)
    speed = upper;
else
    objective = @(v) powerAtVelocity(v,slope*v,energy)/v;
    speed = fminbnd(objective,lower,upper,optimset('TolX',1e-6,'Display','off'));
    candidates = [lower,speed,upper];
    values = arrayfun(objective,candidates);
    [~,which] = min(values);
    speed = candidates(which);
end
speedCache(key) = speed;
end

function speed = optimalVerticalSpeed(direction,lower,upper,energy)
if upper <= lower*(1+1e-9)
    speed = upper;
    return;
end
objective = @(v) powerAtVelocity(0,direction*v,energy)/v;
speed = fminbnd(objective,lower,upper,optimset('TolX',1e-6,'Display','off'));
candidates = [lower,speed,upper];
values = arrayfun(objective,candidates);
[~,which] = min(values);
speed = candidates(which);
end

function power_W = powerAtVelocity(horizontalSpeed,verticalSpeed,energy)
power_W = powerAtDuration(1,horizontalSpeed,verticalSpeed,energy);
end

function power_W = powerAtDuration(time_s,horizontal,vertical,energy)
horizontalSpeed = horizontal/time_s;
verticalSpeed = vertical/time_s;

n = energy.rotorCount;
rho = energy.airDensity;
A = energy.rotorDiscArea;
W = energy.weight;
v0 = energy.inducedVelocity;

inducedPower = (1+energy.inducedCorrection)*W^(3/2)/sqrt(2*n*rho*A);
bladePower = W^(3/2)/sqrt(n*rho*A)* ...
    energy.thrustCoefficient^(-3/2)*energy.profileDrag* ...
    energy.rotorSolidity/8;
hoverPower = bladePower+inducedPower;

inducedRatio = sqrt(1+horizontalSpeed^4/(4*v0^4))- ...
    horizontalSpeed^2/(2*v0^2);
horizontalIncrement = 3/8*sqrt(n)*energy.profileDrag* ...
    sqrt(W*rho*A/energy.thrustCoefficient)*energy.rotorSolidity* ...
    horizontalSpeed^2 + inducedPower*(sqrt(max(inducedRatio,0))-1) + ...
    n/2*energy.flatPlateHorizontal*rho*horizontalSpeed^3;

if abs(verticalSpeed) <= 1e-12
    verticalIncrement = 0;
else
    directionSign = sign(verticalSpeed);
    speedMagnitude = abs(verticalSpeed);
    rootArgument = (1+directionSign*energy.flatPlateVertical/A)* ...
        speedMagnitude^2 + 2*W/(n*rho*A);
    assert(rootArgument >= 0, ...
        'Vertical speed is outside the valid power-model domain.');
    verticalIncrement = 0.5*W*speedMagnitude + ...
        directionSign*n/4*energy.flatPlateVertical*rho*speedMagnitude^3 + ...
        (0.5*W+directionSign*n/4*energy.flatPlateVertical*rho* ...
        speedMagnitude^2)*sqrt(rootArgument);
end

power_W = hoverPower+horizontalIncrement+verticalIncrement;
end
