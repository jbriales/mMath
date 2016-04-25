function hom_x = makehom( x )
% hom_x = makehom( x )
% Returns the homogeneous vector (append 1 at the end of a column vector)

if size(x,2) ~= 1
  error('Use only column vectors with makehom');
end

hom_x = [x;1];