function M = skew(v)
% M = skew(v)
% Creates a skew symmetrix matrix from a vector or list of 3 elements

% This only makes sense for 3-dimensional vectors or lists
assert(numel(v)==3)

if ~iscell(v)
  M = [  0    -v(3)  +v(2)
       +v(3)    0    -v(1)
       -v(2)  +v(1)    0  ];
else
  s  = size(v{1});
  Os = zeros(s);
  M  = [Os -v{3} v{2}; v{3} Os -v{1}; -v{2} v{1} Os];
end

end