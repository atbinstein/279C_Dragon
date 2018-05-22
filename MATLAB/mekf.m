function [xhist, Phist] = mekf(x0, P0, W, V, rN, whist, yhist, dt)

xhist = zeros(7,size(yhist,2));
xhist(:,1) = x0;

Phist = zeros(6,6,size(yhist,2));
Phist(:,:,1) = P0;


for k = 1:(size(yhist,2)-1)
   
    [x_p, A] = prediction(xhist(:,k),whist(:,k),dt);
    yhist(:,k);
    A;
    Phist(:,:,k);
    P_p = A*Phist(:,:,k)*A' + 10*W;
    P_p;
    
    [yp, C] = measurement(x_p(1:4),rN); 
    
    
    %Innovation
    zq = q2phi(qmult(qconj(x_p(1:4)))*yhist(1:4,k+1));
    zr = yhist(5:end,k+1) - yp(5:end);
    z = [zq;zr];
    size(C);
    size(P_p);
    size(V);
    S = C*P_p*C' + V;
%     S = C*P_p*C';

    %Kalman Gain
    L = P_p*C'*S^-1;
    
    %Update
    dx = L*z;
    phi = dx(1:3);
%     dq = [1/2*phi ; 1 - 1/8*phi'*phi];
    dq = phi2q(phi);
    dq = dq/norm(dq);
    
    xhist(1:4,k+1) = qmult(x_p(1:4))*dq; %quaternion update
    
    size(dx);
    xhist(5:7,k+1) = x_p(5:7) + dx(4:6); %bias update
    P_p;
    L*C;
    Phist(:,:,k+1) = (eye(6) - L*C)*P_p*(eye(6) - L*C)' + L*V*L'; %covariance update
    x_p(5:7);

end


end

