function [xn,A] = prediction(xk,w,dt)

q = xk(1:4);
b = xk(5:7);
theta = norm(w-b)*dt;
r = (w-b)/norm(w-b);
xn(1:4) = qmult(q)*[r*sin(theta/2);cos(theta/2)];
xn(5:7) = b;
xn = xn';
R = expm(hat(-w+b)*dt);

A = [R -dt*eye(3) ; zeros(3) eye(3)];

end

