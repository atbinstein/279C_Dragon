function A = hat(e)
 
A = zeros(3,3);
A(1,2) = -e(3);
A(2,1) = e(3);
A(1,3) = e(2);
A(3,1) = -e(2);
A(2,3) = -e(1);
A(3,2) = e(1);

end