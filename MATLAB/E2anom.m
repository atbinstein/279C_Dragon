function anom = E2anom(E, e)
% E2anom Computes true anomaly given eccentric anomaly and eccentricity
%
% Inputs:
%       E - eccentric anomaly [rad]
%      e - eccentricity of orbit
%
% Outputs:
%   anom - mean anomaly [rad]

anom = acos((cos(E) - e)/(1 - e*cos(E)));

% Make sure E sits in the correct semi-plane
if E > pi
    anom = 2*pi - anom;
end

end