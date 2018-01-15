function out = isrotation(X,tol)
% is = isrotation(X)
% is = isrotation(X,tol)
% Check if the input is a rotation matrix

if nargin < 2
  tol = 1e-4;
end

if norm(X'*X-eye(3),'fro') > tol || norm(det(X)-1) > tol
  out = false;
else
  out = true;
end