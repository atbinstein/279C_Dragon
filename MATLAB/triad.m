function [q,R] = triad(rN1,rN2,rB1,rB2)

tn1 = rN1;
tn2 = cross(rN1,rN2)/norm(cross(rN1,rN2));

Tn = [tn1 tn2 cross(tn1,tn2)/norm(cross(tn1,tn2))];

tb1 = rB1;
tb2 = cross(rB1,rB2)/norm(cross(rB1,rB2));

Tb = [tb1 tb2 cross(tb1,tb2)/norm(cross(tb1,tb2))];

R = Tn*Tb';
 
q = Qtoq(R);

end