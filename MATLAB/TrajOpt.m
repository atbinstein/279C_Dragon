%% Attitude Trajectory Optimization
% Arthur Binstein and Marco Hinojosa
load EigenStates.mat

global J N; 
J = diag(rand(1,3)*100);
N = size(qOpt,2);

obj = @(x) sum(x(1:N));

x0 = [1*ones(N,1); reshape(qOpt,4*N,1); reshape(omOpt,3*N,1); reshape(tauOpt,3*N,1)];

A = [];
b = [];
Aeq = [];
beq = [];

lb = [zeros(N,1); -Inf*ones(4*N,1); -Inf*ones(3*N,1); -Inf*ones(3*N,1)];
ub = [Inf*ones(N,1); Inf*ones(4*N,1); Inf*ones(3*N,1); Inf*ones(3*N,1)];

options = optimoptions(@fmincon, 'TolFun', .000000001, 'MaxIter', 10000, ...
                       'MaxFunEvals', 100000, 'Algorithm', 'sqp');

optimal = fmincon(obj,x0,A,b,Aeq,beq,lb,ub,@dynamics, options);


function [c, ceq] = dynamics(x)
global N J;
qd = [1 2 3 4]';
qd = qd/norm(qd);
c = [];
tf = sum(x(1:N));
q = reshape(x(1+N:5*N),4,N);
w = reshape(x(1+5*N:8*N),3,N);
u = reshape(x(1+8*N:end),3,N);

dt = tf/N;
ceq = [q(:,1) - [0 0 0 1]'; w(:,1); norm(q(:,1)) - 1];
for ii = 1:N-1
    x_i = [q(:,ii); w(:,ii)];
    x_n = [q(:,ii+1); w(:,ii+1)];
    xdot_i = [.5*qmult(q(:,ii))*[w(:,ii); 0]; 
               -inv(J)*(cross(w(:,ii),J*w(:,ii)) - u(:,ii))];
    xdot_n = [.5*qmult(q(:,ii+1))*[w(:,ii+1); 0]; 
               -inv(J)*(cross(w(:,ii+1),J*w(:,ii+1)) - u(:,ii+1))];
           
    xend = x_i + dt/2*(xdot_i + xdot_n);
    ceq = [ceq; x_n - xend; norm(q(:,ii)) - 1; x(ii+1) - x(ii)];
end
ceq = [ceq; q(:,end) - qd; w(:,end); norm(q(:,end)) - 1]; 
end