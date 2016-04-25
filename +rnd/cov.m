function A = cov(n)
% Generate random nxn symmetric positive semidefinite matrix
%   A = cov(n)
% This is useful for getting covariance matrices

[U,S,~] = svd(randn(n));
A = U*S*U';