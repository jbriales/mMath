A collection of utilities related to linear algebra.

## Common *matrix* vector spaces

These are vector spaces whose elements may be represented by matrices
(similarly to the concept of *matrix manifolds*).

We consider the following common vector spaces
(note the associated *Tag* between parentheses):
- `Re^{m x n}` (*Re*): Rectangular matrices
- `Sym(d)` (*Sym*): Symmetric dxd matrices
- `Skew(d)` (*Skew*): Skew-symmetric or anti-symmetric dxd matrices

The general operations or methods available for each space are:
- `v = vec(M)`: vectorize an element M in the space S.
  In ReSpace this is just a vector concatenating all elements in the matrix column-wise.
  In general, this returns a minimal representation of the element in terms of a *vector of coordinates*
  corresponding to the basis elements.
  The returned vector is a *column* vector.
- `M = mat(v,dim)`: un-vectorize a vector of coordinates into the corresponding element in the space S.
  This is the inverse operation of `vec(M)`.
- `dimVec`: inherent dimension of the vector space (dimension of its basis)
- `dimMat`: dimension of the matrix elements (ambient space)
- `basis`: returns a cell array with canonical basis elements

For each space and operation, we use the camelCase syntax `<methodTag><SpaceTag>`.

Common shortcuts:

Back-compatibility commands:

Other utilities:
- `symPart`
- `skewPart`

X = mat(x,n)


