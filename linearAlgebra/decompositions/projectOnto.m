function Xproj = projectOnto(U,X)
% Xproj = projectOnto(U,X)
% Project column vectors in X onto the subspace spanned by U columns
Xproj = U*(U'*X);