function e = canvec( n, k )

assert( k >= 1 && k <= n );

In = eye(n);
e = In(:,k);

end