# Mmath toolbox
Collection of Matlab math-related functions for improved readability and compactness.

- `R = project_rotation( M )`
  Find rotation `R=U*S*V'` closest to `M=U*D*V'` minimizing the chordal distance.

- `[R,t] = project_rotation( M, fmin )`
  Find rotation `R=U*S*V'` that minimizes `fmin(R)`.

- `Z = null(A,tol)`
  The same function as Matlab but with user-definable tolerance.

- `Q = orth(A,tol)`
  The same function as Matlab but with user-definable tolerance.
