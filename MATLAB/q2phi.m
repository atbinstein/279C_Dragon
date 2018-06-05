function phi = q2phi(q)

% r = q(1:3)/norm(q(1:3));
% theta = 2*acos(q(4));
% 
% phi = r*theta;

Q = q2Q(q);

phi = unhat(logm(Q));

end