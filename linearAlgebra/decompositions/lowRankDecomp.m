function [U,sv] = lowRankDecomp(M,rnk)
% U = lowRankDecomp(M)
% Returns the range orthogonal basis so that U*U'=M.
%
% U = lowRankDecomp(M,rnk)
% Keep only the rnk largest singular vectors.

if nargin == 2
  [U,S,~] = svd(M);
  U = U(:,1:rnk);
else % nargin == 1
  U = orth(M);
  error('Set Singular Values')
end
sv = diag(S);

end
