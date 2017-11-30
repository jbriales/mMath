function x = vecs(X)
% x = vecs(X)
% Symmetric (half) vectorization for symmetric matrices
% This vectorization is scaled so that tr(ST)=vecs(S)'*vecs(T)
% Taken from "Schäcke, K. (2013). On the kronecker product."
% 
% See also: SymSpace(...).vecTr

s = size(X);

% Scale the off-diagonal lower triangular block with sqrt(2)
X( ~triu(true(s)) ) = sqrt(2) * X( ~triu(true(s)) );

% Extract lower triangular block to vector
x = X( tril(true(s)) );

end