function rB = noisyMagno(rN,q,E,mu)

phi = chol(E)*randn(3,1); %Random error as axis-angle

Qbn = q2Q(q);
Qnb = Qbn';
rB = [];
for ii = 1:size(rN,2)
    Qnb*rN(:,ii);
    
end



end