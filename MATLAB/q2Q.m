function Q = q2Q(q)
v = q(1:3);
s = q(4);

Q = eye(3) + 2*hat(v)*(s*eye(3) + hat(v));

end