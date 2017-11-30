function S = symPart(M)
% sym(M) = symPart(M)
% Extract the symmetric part of a matrix,
%   A = 0.5*(A+A')  +  0.5*(A-A'),
% according to the usual direct sum relation
%   Re^{dxd} = Sym(d) \oplus Skew(d).
% 
% See also: skewPart, symmetrize
S = 0.5*(M+M');
end