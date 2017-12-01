function M = blkskew(m1,m2,m3)
% M = blkskew(m1,m2,m3)
% Creates a skew symmetrix matrix from 3 blocks as
% 
% See also SKEW.

error('TODO: Implement this properly')
% This is the original code in skew for cells
s  = size(v{1});
Os = zeros(s);
M  = [ Os   -v{3}  v{2};
  v{3}  Os   -v{1};
  -v{2}  v{1}   Os];