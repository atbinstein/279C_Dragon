function M = E2M(E, e)
% E2M Computes mean anomaly given eccentric anomaly and eccentricity
%
% Inputs:
%     E - eccentric anomaly [rad]
%     e - eccentricity of orbit
%
% Outputs:
%     M - mean anomaly [rad]

% Kepler's equation
M = E - e*sin(E);

end