function e = canvec( n, k )
% e = canvec( n, k )
% Returns the canonical col-vector of dimension n with 1 at k-th position.

assert( k >= 1 && k <= n );

% In = eye(n);
% e = In(:,k);
e = full( sparse(k,1, 1, n,1, 1) );

end