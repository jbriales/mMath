function A = cov(n)
% Generate random nxn symmetric positive semidefinite matrix
%   A = cov(n)
% This is useful for getting covariance matrices

[U,S,~] = svd(randn(n));
A = U*S*U';

% To avoid numerical issues, the result is symmetrized!
% Otherwise epsilon magnitude errors appear in norm(A-A')
A = symmetrize(A);