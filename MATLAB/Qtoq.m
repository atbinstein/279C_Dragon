function q = Qtoq(Q)

phi = unhat(logm(Q));

theta = norm(phi);

r = phi/theta;

q = [r*sin(theta/2) ; cos(theta/2)];

end
