classdef LinMatSpace
  %LinMatSpace
  %   This is the abstract base class for vector spaces of interest
  %   whose elements may be represented by matrices (akin to matrix
  %   manifolds).
  %
  % Common properties:
  %   - dimVec: dimension of vector space (set of coordinates)
  %   - dimMat: dimension of ambient space (matrix or array)
  %   - basis: cell array of basis elements
  % Common methods:
  %   - v = vec(M)
  %   - M = mat(v)
  %   
  % See also 
  
  properties( SetAccess = protected )
    vdim % dimension of the vector space (dim of vector basis)
    Mdim % dimension of the matrix representation (ambient space)
  end
  properties( Dependent )
    dimVec
    dimMat
  end
  methods
    function d = get.dimVec(this), d = this.vdim; end
    function D = get.dimMat(this), D = this.Mdim; end
  end
  properties( Dependent, Abstract )
    basis
  end
  
  methods 
    function this = LinMatSpace( dims )
      this.Mdim = dims;
    end
  end
  
  methods (Abstract)
    v = vec(this, M)
    M = mat(this, v)
%     d = dimVec( )
%     D = dimMat( )
%     c_M = basis( )
  end
  
end

