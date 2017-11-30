function X = vec_inv(x)
% DEPRECATED. Use M = mat(v) instead.
if nargin == 1
  dim = sqrt(numel(x));
  assert( dim-fix(dim) == 0 ); % Integer value
  X = reshape(x,dim,dim);
end

end