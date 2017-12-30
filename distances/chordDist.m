function d = chordDist(X,Y)
% d = chordDist(X,Y)
% Compute chordal distance between two arrays, given by Frobenius norm 
% of difference.
% 
% d = chordDist(R1,R2)
% For rotations, this distance is related to the angular distance by
%   chordDist = 2\sqrt(2)\sin(theta/2)
% 
% See also angDist.

d = norm(X-Y,'fro');

end