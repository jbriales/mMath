function A = stackArr(varargin)
% [A] = stackArr(A1,...,An)
% Stack homogeneous set of array elements into next empty dimension.
% Note our convention is to use column vectors, so stack occurs
% - in DIM=2 for sets of column vectors,
% - in DIM=3 for sets of matrices,
% - so on for higher dimensions if necessary.
% 
% See also unstackArr.

% Get next non-used dimension
A1 = varargin{1};
if isvector(A1)
  assert(size(A1,2)==1)
  stackDim = 2;
else
  stackDim = ndims()+1;
end
% Concatenate all elements along this dimension
A = cat(stackDim,varargin{:});