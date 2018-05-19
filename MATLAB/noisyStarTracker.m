function y = noisyStarTracker(q,phi)

qq = phi2q(phi);
y = qmult(q)*qq;

end
