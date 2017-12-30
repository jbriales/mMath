function [dmin,Xmin] = closestDist(Xref,c_X, distFun)
% [dmin,Xmin] = closestDist(Xref,c_X, distFun)
% Find the element in the cell array c_X closest to Xref, according to
% the distance function distFun (default: chordDist).
% Return minimum distance and corresponding element.

if nargin < 3
  % Default: chordal distance of two arrays (Frobenius norm)
  distFun = @chordDist;
end

assert(iscell(c_X))
N = numel(c_X);
d = NaN(1,N);
for k=1:N
  d(k) = distFun(Xref,c_X{k});
end
[dmin,kmin] = min(d);
Xmin = c_X{kmin};

end