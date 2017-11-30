function x = vech(X)
% x = vech(X)
% 
% Half vectorization for symmetric matrices
% Take only elements in the lower triangular block of the matrix
% 
% See also: SymSpace(...).vec method

s = size(X);
x = X( tril(true(s)) );

end