function [dmin,Xmin] = closestDist(Xref,l_X, distFun)
% [dmin,Xmin] = closestDist(Xref,l_X, distFun)
% Find the element in a list l_X closest to Xref, according to
% the distance function distFun (default: chordDist).
% Return minimum distance and corresponding element.
% The list can be:
% - A cell array c_X of any element
% - An array a_X of objects

if nargin < 3
  % Default: chordal distance of two arrays (Frobenius norm)
  distFun = @chordDist;
end

assert(~isempty(l_X),'List should not be empty')
if isobject(l_X)
  % Force into cell array to unify code
  l_X = num2cell(l_X);
end
assert(iscell(l_X))
c_X = l_X;
N = numel(c_X);
d = NaN(1,N);
for k=1:N
  d(k) = distFun(Xref,c_X{k});
end
[dmin,kmin] = min(d);
Xmin = c_X{kmin};

end