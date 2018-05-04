function e = unhat(A)

e = zeros(3,1);

e(1) = -A(2,3);
e(2) = A(1,3);
e(3) = -A(1,2);

end