function c_A = unstackArr(A)
% c_A = unstackArr(DIM,[A])
% Divide stacked layers of A and return cell array.
% Note our convention is to use column vectors, so stack occurs
% - in DIM=2 for sets of column vectors,
% - in DIM=3 for sets of matrices,
% - so on for higher dimensions if necessary.
% 
% See also stackArr.

% Get maximum dimension of array (the one used for stacking)
stackDim = ndims(A);
% Get dimensions of each singular element stacked
elemDims = 1:stackDim-1;

% Divide stacked array into cell of elements
c_A = num2cell(A,elemDims);
% The direction of num2cell leaves singleton dimensions, squeeze!
c_A = squeeze(c_A);
% Prefer row cell arrays for sets
c_A = transpose( c_A(:) );