function d = antiDist(q1,q2)
% d = antiDist(q1,q2)
% Compute a form of chordal distance which is antipodally symmetric
% on the n-sphere:
%   d(X,Y) = 1 - <X,Y>^2
% Thanks to the square, the distance is unaffected by sign of X or Y.
% For normalized vectors, this distance should fulfill 0 < d < 1.
% 
% This distance is also connected to the Bingham distribution,
% and for quaternion vectors it should fulfill
%   d(q1,q2) = 1 - <q1,q2>^2 = (1-cos(theta))/2
% Ref: https://math.stackexchange.com/a/90098/344323
% 
% See also chordDist.

d = 1 - dot(q1,q2)^2;
% d = 1 - trace(X'*Y)^2;

end