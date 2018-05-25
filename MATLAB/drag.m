function [D,T] = drag(c,n,A,V,r,q)
r = r*1000;
V = V*1000;
R_e = 6378*1000; %m
h0 = 0;
h = norm(r) - R_e;
p0 = 1.225; %kg/m^3
H = 15*1000; %m
% rho = p0*exp(-(h - h0)/H);
rho = 2.72*10^-12;
nQb = q2Q(q);
omE = .000072921158553*[0 0 1]'; %rad/s
Vn = (V + cross(omE,r)); %m/s
Vb = nQb'*Vn; 
Cd = 2.2;

F = zeros(size(c));
D = zeros(3,1);
Tau = zeros(size(c));
T = zeros(3,1);
for ii = 1:size(c,2)
    F(:,ii) = .5*rho*Cd*norm(Vb)*Vb*A(ii)* ...
        max(dot(Vb,n(:,ii))/norm(Vb),0);
    D = D + F(:,ii);
    Tau(:,ii) = cross(c(:,ii),F(:,ii));
    T = T + Tau(:,ii);
end