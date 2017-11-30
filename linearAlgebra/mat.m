function M = mat(v,dims)
% M = mat(v) or M = mat(v,dims) (the 2nd argument is optional)
%   Given a vector of length d1 ... dn, this produces the n x n matrix
%   Y such that x = vec(Y).  In other words, x contains the columns of the
%   matrix Y, stacked below each other.
% 
% See also vec.

% Adapted from Sedumi

if nargin < 2
  % Assume square matrix and try
  d = floor(sqrt(length(v)));
  if (d*d) ~= length(v)
    error('Argument X has to be a square matrix')
  end
  dims = [d d];
end

M = reshape(v,dims);