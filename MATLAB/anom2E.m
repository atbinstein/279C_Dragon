function E = anom2E(anom, e)
% anom2E Computes eccentric anomaly given true anomaly and eccentricity
%
% Inputs:
%   anom - mean anomaly [rad]
%      e - eccentricity of orbit
%
% Outputs:
%      E - eccentric anomaly [rad]

E = acos((e + cos(anom))/(1 + e*cos(anom)));

% Make sure E sits in the correct semi-plane
if anom > pi
    E = 2*pi - E;
end

end
