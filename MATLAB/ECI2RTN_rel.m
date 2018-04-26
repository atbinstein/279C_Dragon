
function [T] = ECI2RTN_rel(r,v)

h = cross(r,v);

Rhat = r/norm(r);
Nhat = h/norm(h);
That = cross(Nhat,Rhat)/norm(cross(Nhat,Rhat));

T = [Rhat,...
   That,...
   Nhat];

end

