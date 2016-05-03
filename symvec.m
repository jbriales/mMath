function x = symvec(var,n,varargin)
% x = symvec(var,n)
% Create a symbolic vector of basename var with n elements
% 
% symmat(___, set) sets the assumption that the variables % belong to a set.
% Here, set can be 'real', 'positive', 'integer', or 'rational'.

x = sym([var,'_%d'],[n 1],varargin{:});