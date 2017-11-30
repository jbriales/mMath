% Script to test spaces

% Check basis for each space
d = 3;
c_spaces = {Re1DSpace(d), Re2DSpace(d), SymSpace(d), SkewSpace(d)};
  
for k = 1:numel(c_spaces)
  S = c_spaces{k};
  
  c_E = S.basis;
  c_E{:}
end
