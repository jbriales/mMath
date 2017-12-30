function coeff = convexCombDecomp( Z0, c_Z )
% coeff = convexCombDecomp( Z0, c_Z )
% Given a target element Z0 and a set of basis elements c_Z,
% find vector of coefficients coeff so that the convex combination fulfills
%   Z0 = sum_k( ak*Zk )

% We treat this as a system of linear equations,
% then check constraints on coefficients a-posteriori.

% Get number of coefficients
% d = numel(d_Z);

% Build data elements for linear system
% Vectorize the matrix vectors
%  Z0 = sum_k( ak*Zk ) -> vec(Z0) = sum_k( vec(Zk)* ak ) = [vec(Zk)]*a
M = [];
for k=1:numel(c_Z)
  Zk = c_Z{k};
  M(:,k) = vec(Zk);
end
% Get independent term
b = vec(Z0);
% Solve convex coefficients
coeff = M \ b;

% Check constraints on coefficients
assert( all(coeff>0) )
assert( sum(coeff)==1 )
