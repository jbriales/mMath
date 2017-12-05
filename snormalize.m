function X_inS = snormalize( X )
% x_inS = snormalize( x )
% Returns the normalized vector so that norm(x_inS)=1
% Works with column vectors
%
% See also hnormalize.

s = size(X);
vec1 = ones(s(1),1);
X_inS = X ./ (vec1*sqrt(sum(X.^2,1)));