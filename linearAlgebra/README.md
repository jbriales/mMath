A collection of utilities related to linear algebra.

## Common *matrix* vector spaces

These are vector spaces whose elements may be represented by matrices
(similarly to the concept of *matrix manifolds*).
We develop a common abstract base class `LinMatSpace` that unifies
the common interface and properties for all particular instances.

Properties:
- `dimVec`: inherent dimension of the vector space (dimension of its basis)
- `dimMat`: dimension of the matrix elements (ambient space)
- `basis`: returns a cell array with canonical basis elements

Methods:
- `v = vec(M)`: vectorize an element M in the space S.
  In ReSpace this is just a vector concatenating all elements in the matrix column-wise.
  In general, this returns a minimal representation of the element in terms of a *vector of coordinates*
  corresponding to the basis elements.
  The returned vector is a *column* vector.
- `M = mat(v,dim)`: un-vectorize a vector of coordinates into the corresponding element in the space S.
  This is the inverse operation of `vec(M)`.

We consider the following common vector spaces and corresponding classes
(note the associated *Tag* between parentheses):
- `Re^{nrows x 1}` (*Re1D*): Column vector space
- `Re^{nrows x ncols}` (*Re2D*): Rectangular matrices
- `Re^{d1 x ... x dN}` (*ReND*): Multi-dimensional array
- `Sym(d)` (*Sym*): Symmetric dxd matrices
- `Skew(d)` (*Skew*): Skew-symmetric or anti-symmetric dxd matrices

For each space and operation, we use the camelCase syntax `<methodTag><SpaceTag>`.

Common shortcuts (kept for compatibility and convenience):
- `canvec`
- `canvec3`
- `canvec33`

Back-compatibility commands:
- `vecs`, `vech`, `avecs`, `avech`, `vec_inv` functions.
TODO:
- Remove these confusing functions where used, use classes instead.

Other utilities:
- `symPart`
- `skewPart`








