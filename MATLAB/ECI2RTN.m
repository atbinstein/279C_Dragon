function [r_rtn,v_rtn] = ECI2RTN(r_eci, v_eci)

Rhat = r_eci/norm(r_eci);
Nhat = cross(r_eci,v_eci)/norm(cross(r_eci,v_eci));
That = cross(Nhat,Rhat)/norm(cross(Nhat,Rhat));


T = [Rhat,...
   That,...
   Nhat];

r_rtn = T*r_eci;
v_rtn = T*v_eci;
end

