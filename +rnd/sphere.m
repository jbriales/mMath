function x = sphere(d,n,R)
% Points in the surface of a sphere in R^d
% x = sphere(d,n,R)
% Parameters:
%   d - dimension of the ambient space: x \in Re^d
%   n - number of points to generate
%   R - radius of the sphere: ||x||_2 == R

if nargin < 3
  R = 1;
end
if nargin < 2
  n = 1;
end
if nargin < 1
  d = 3;
end

% points in the R-(d-1)-sphere
x = R*snormalize(randn(d,n));

end
