function q = phi2q(phi)

if norm(phi) > 10*pi/180
    theta = norm(phi);
    r = phi/theta;

    q(1:3) = r.*sin(theta/2);
    q(4) = cos(theta/2);
    q = q';
else
    q = [phi./2 ; 1 - 1/8*phi'*phi];
    q = q/norm(q);
end
end