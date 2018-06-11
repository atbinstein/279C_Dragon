%% Attitude Trajectory Optimization
% Arthur Binstein and Marco Hinojosa
close all; clear all;
load EigenStates.mat

% J = diag([20000 50000 60000]);
N = size(qOpt,2);

x0 = [qOpt; omOpt; tauOpt; dt*ones(1,N)];
x0 = x0(:);
x0(end-3:end) = [];

Aeq = zeros(N-1,11*N-4);
c = 1;
for ii = 1:N-2
    Aeq(ii,11*ii) = 1;
    Aeq(ii,11*ii+11) = -1;
end
Aeq = [eye(7), zeros(7,11*(N-1)); zeros(7,11*(N-1)) eye(7); Aeq];
% Aeq = sparse(Aeq);

objGrad = zeros(1,11*N-4);
for ii = 1:length(objGrad)
    if mod(ii,11) == 0
        objGrad(ii) = 1;
    end
end

beq = zeros(N-1,1);
beq = [x0(1:7); qOpt(:,end); omOpt(:,end) ; beq];
A = [];
b = [];

lb = [-Inf*ones(4,N); -100*ones(3,N); -1500*ones(3,N); dt/10*ones(1,N)];
ub = [Inf*ones(4,N); 100*ones(3,N); 1500*ones(3,N); dt*1.5*ones(1,N)];
lb = lb(:);
ub = ub(:);
lb(end-3:end) = [];
ub(end-3:end) = [];

options = optimoptions(@fmincon, 'TolFun', .001, 'MaxIter', 10000, ...
    'MaxFunEvals', 100000, 'Algorithm','sqp',...
    'Display','iter','SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient',true,'ConstraintTolerance',.001);
% options = optimoptions(@fmincon, 'TolFun', .0001, 'MaxIter', 10000, ...
%     'MaxFunEvals', 100000, 'Algorithm', 'sqp', ...
%     'Display','iter','SpecifyObjectiveGradient',true);


optimal = fmincon(@(x)fun(x,objGrad),x0,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,N,J,q0,qd), options);

function [obj, objGrad] = fun(x,objGrad)

obj = objGrad*x;
objGrad = objGrad;

end

function [c, ceq, dc, dceq] = nonlcon(x,N,J,q0,qd)
% function [c, ceq] = nonlcon(x,N,J,q0,qd)
c = [];
dc = [];

xmod = reshape(x(1:11*(N-1)),11,N-1);
q = [xmod(1:4,:) x(end-6:end-3)];
w = [xmod(5:7,:) x(end-2:end)];
u = xmod(8:10,:);
dt = xmod(11,:);

% [fx, dfx, du] = dynamics(q(:,1),w(:,1),u(:,1),J);
% A = fillrows(fx, dfx, du, dt(1),q(:,1),0,length(x));
% dceq = A;
ceq = zeros(7*(N-1)+N-2,1);
dceq = zeros(7*(N-1)+N-2,11*N-4);
% size(dceq)
for ii = 1:N-1
    ii;
    qm = .5*(q(:,ii) + q(:,ii+1));
    wm = .5*(w(:,ii) + w(:,ii+1));
    [fx, dfx, du] = dynamics(qm,wm,u(:,ii),J);
    ceq(1+8*(ii-1):7+8*(ii-1),1) = [q(:,ii); w(:,ii)] + fx*dt(ii) - [q(:,ii+1); w(:,ii+1)];
    if ii ~= (N-1)
        ceq(8*ii,1) = q(:,ii)'*q(:,ii) - 1;
    end
    size(ceq);
    A = [];
    A = fillrows(fx, dfx, du, dt(ii), q(:,ii), ii, length(x));
    size(A);
    size(ceq);
    size(dceq);
    if ii ~= (N-1)
        dceq(1+8*(ii-1):8+8*(ii-1),:) =  A;
    else
        dceq(1+8*(ii-1):7+8*(ii-1),:) =  A;
    end
end
dceq = dceq';
end

function [fx,dfx,du] = dynamics(q,w,u,J)
% function [fx] = dynamics(q,w,u,J)

fx(1:4,1) = .5*qmult(q)*[w;0];
fx(5:7,1) = -inv(J)*(cross(w,J*w) - u);

dfx = [.5*qmultR([w;0]) .5*qmult(q)*[eye(3); zeros(1,3)];
       zeros(3,4)       -inv(J)*(hat(w)*J - hat(J*w))];

du = [zeros(4,3); inv(J)];

end

function A = fillrows(fx, dfx, du, dt, q, n, N)
m = n-1;

if n ~= (N+4)/11 - 1
    A = zeros(8,N);
    A(1:7,m*11+1:m*11+7) = eye(7) + dfx*.5*eye(7)*dt;
    A(1:7,m*11+8:m*11+10) = dt*du;
    A(1:7,m*11+11) = fx;
    A(1:7,m*11+12:m*11+18) = dfx*.5*eye(7)*dt - eye(7);
    A(8,m*11+1:m*11+4) = 2*q';  
else
    A = zeros(7,N);
    A(1:7,m*11+1:m*11+7) = eye(7) + dfx*.5*eye(7)*dt;
    A(1:7,m*11+8:m*11+10) = dt*du;
    A(1:7,m*11+11) = fx;
    A(1:7,m*11+12:m*11+18) = dfx*.5*eye(7)*dt - eye(7);
end

end

