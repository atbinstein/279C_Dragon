function qhat = qmultR(q)

v = q(1:3);
s = q(4);

qhat = [s*eye(3) - hat(v)  v ; -v' s]; 

end