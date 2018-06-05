function [y,beta] = addNoise(x,V,flag,rN,beta)
% HSTmeasure: Creates sensor measurements by adding noise to true states.
%
% Inputs:
%      x - Ttime History of True States
%      V - Sensor Noise Covariance
%      flag - Sensor Type
%
% Outputs:
%   y - Noisy Sensor Measurements
%   beta - Bias History (for Gyros)

if isempty(rN)
else
    rN1 = rN(:,1);
    rN2 = rN(:,2);
    rN3 = rN(:,3);
    beta = [];
end

if strcmp(flag,'Magno') % Magnetometers Body Vectors
    q_stat = x(4:7,1)';
    Q = q2Q(q_stat);

    rB1(:) = (Q'*rN1) + chol(V)*randn(3,1);
    rB1(:) = rB1(:)/norm(rB1(:));

    y(:) = rB1(:);

elseif strcmp(flag,'StarTracker') % Star Tracker (Quaternion)
    phi = chol(V)*randn(3,1); % Axis-Angle Noise
    dq = phi2q(phi); % Convert to Quaternion Noise
    dq = dq/norm(dq); % Re-Normalize
    y(:) = qmult(x(4:7,1))*dq; % "Add" Noise by Multiplying to State Quaternion

elseif strcmp(flag,'GPS')
    q_stat = x(4:7,1)';
    Q = q2Q(q_stat);

    rB2(:) = (Q'*rN2) + chol(V)*randn(3,1);
    rB2(:) = rB2(:)/norm(rB2(:));
    rB3(:) = (Q'*rN3) + chol(V)*randn(3,1);
    rB3(:) = rB3(:)/norm(rB3(:));
    y(:,1) = [rB2' ; rB3'];
    
elseif strcmp(flag,'Gyro') % Gyroscope (Vector)
    beta(:,1) = [0;0;0];
    W_rnd = V(1:3,1:3); % Rate Noise Density Covariance
    W_arw = V(4:6,4:6); % Angle Random Walk Covariance
    om_true = x(1:3,1);
    beta = beta + chol(W_arw)*randn(3,1);
    om_noise(:,1) = om_true + beta + chol(W_rnd)*randn(3,1);
    y(:,1) = om_noise;
end
end