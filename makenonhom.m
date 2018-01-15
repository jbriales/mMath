function nonhom_x = makenonhom( x )
% nonhom_x = makenonhom( x )
% Returns the non-homogeneous vector: normalize dividing by the last
% element in a column vector and drop it

if any( abs(x(end,:)) < 1e-6 )
  warning('Hom component is zero');
end

if isvector(x)
  nonhom_x = x(1:end-1) / x(end);
else
  n = size(x,1) - 1;
  nonhom_x = x(1:n,:) ./ repmat(x(end,:),n,1);
end