function q = qconj(qq)

v = qq(1:3);
v = -v;

q = [v; qq(4)];

end
