function [y,C] = measurement(q,rN)

QBN = q2Q(q);

rB = [];
for ii = 1:size(rN,2)
    rB = [rB; QBN'*rN(:,ii)];
end
rB;

C = [eye(3) zeros(3)];
for ii = 1:size(rN,2)
    C = [C ; hat(rB(ii*3-2:ii*3)) zeros(3)];
end

y = [q; rB(:)];

end