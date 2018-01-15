function [R,f] = projectSO3( M, fmin )
% R = projectSO3( M )
% 
% Takes a general linear matrix M and projects it to
% the rotation matrix in SO(3) minimizing the chordal distance.
% This is equivalent to find the closest matrix under Frobenius norm
% (see "Hartley et al. Rotation averaging" in IJCV).
% 
% The problem is solved using the SVD decomposition of M:
%   M = U*S*V'
% If det(U*V')>0, R = U*V'
% BUT if det(U*V')<0, R = U*diag([+1 +1 -1])*V'

% Project to valid rotation
[U,D,V] = svd(M);

if det(U*V')>0
  R = U*V';
else
  R = U*diag([1 1 -1])*V';
end

end
