function [d,res] = distanceTo(U,X)
% d = distanceTo(U,X)
% [d,res] = distanceTo(U,X)
% Measures the distance of a point to a subspace.

res = X - projectOnto(U,X);
d = norm(res,'fro');