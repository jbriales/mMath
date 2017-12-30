function d = eucDist(X,Y)
% d = eucDist(X,Y)
% Compute Euclidean distance between two arrays, given by 2-norm
% of difference.
% 
% See also chordDist.

d = norm(X-Y,2);

end