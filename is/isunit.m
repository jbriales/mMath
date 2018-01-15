function out = isunit(X,tol)
% is = isrotation(v)
% is = isrotation(v,tol)
% Check if the input is a unit vector.
% It also works for column vectors stacked horizontally.

if nargin < 2
  tol = 1e-6;
end

s = size(X,2);
out = ( abs(sqrt(sum(X.^2,1)) - ones(1,s)) < tol );

end