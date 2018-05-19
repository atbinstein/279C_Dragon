function y = noisyMagno(q,E,mu)

phi = chol(E)*randn(3,1);

qq = phi2q(phi);

y = qmult(q)*qq;

end