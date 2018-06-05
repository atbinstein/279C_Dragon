function [x, P] = mekfPID(xk, Pk, W, V, rN, wk, ykk, dt)
xk;
wk;
dt;
[x_p, A] = prediction(xk,wk,dt);
P_p = A*Pk*A' + 10*W;

[yp, C] = measurement(x_p(1:4),rN);


%Innovation
x_p(1:4);
ykk(1:4);
zq = q2phi(qmult(qconj(x_p(1:4)))*ykk(1:4));
zr = ykk(5:end) - yp(5:end);
z = [zq;zr];
size(C);
size(P_p);
size(V);
S = C*P_p*C' + V;
%     S = C*P_p*C';

%Kalman Gain
L = P_p*C'*S^-1;

%Update
dx = L*z;
phi = dx(1:3);
%     dq = [1/2*phi ; 1 - 1/8*phi'*phi];
dq = phi2q(phi);
dq = dq/norm(dq);

x(1:4,1) = qmult(x_p(1:4))*dq; %quaternion update

size(dx);
x(5:7,1) = x_p(5:7) + dx(4:6); %bias update
P_p;
L*C;
P = (eye(6) - L*C)*P_p*(eye(6) - L*C)' + L*V*L'; %covariance update
x_p(5:7);

end

