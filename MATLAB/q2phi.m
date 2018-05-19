function phi = q2phi(q)

Q = q2Q(q);

phi = unhat(logm(Q));

end