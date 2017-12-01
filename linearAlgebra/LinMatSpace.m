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
  % Vectorization and its inverse methods:
  %   - v = vec(M)
  %   - M = mat(v)
  % Linear index and subscript correspondences for basis elements:
  %   - [i1,...,iN] = ind2sub(k)
  %   - k = sub2ind(i1,...,iN)
  %   
  % See also 
  
  properties( GetAccess = protected )
    % user cannot see these properties (internal use only)
    vdim % dimension of the vector space (dim of vector basis)
    Mdim % dimension of the matrix representation (ambient space)
    m_ij % matrix of subscripts (each row the k-th subscripts)
    m_k  % matrix of linear indeces (each ij element has k index)
  end
  properties( SetAccess = protected )
    % user can read but not write these properties
    basis
  end
  properties( Dependent )
    dimVec
    dimMat
  end
  methods
    function d = get.dimVec(this), d = this.vdim; end
    function D = get.dimMat(this), D = this.Mdim; end
  end
    
  methods 
    function this = LinMatSpace( dims )
      % For abstract inheriting classes
      if nargin == 0
        return
      end
      
      this.Mdim = dims;
      % NOTE: The constructor should take care of setting properties
      % - basis
      % - m_ij
      % - m_k
      % Initialize
      this.basis = cell(1,this.dimVec);
      this.m_ij = zeros( this.dimVec, numel(dims) );
      this.m_k  = zeros( dims );
    end
    
    % Conversion methods between linear *ind*ex and array *sub*scripts
    function [varargout] = ind2sub(this, k)
      subscripts = this.m_ij(k,:);
      c_ij = num2cell( subscripts );
      varargout = c_ij;
    end
    function k = sub2ind(this, varargin)
      k = this.m_k(varargin{:});
    end
  end
  
  methods (Abstract)
    % Vectorization operation: get vector coordinates for element
    v = vec(this, M)
    % Inverse of vectorization: build mat element from vector coordinates
    M = mat(this, v)
%     d = dimVec( )
%     D = dimMat( )
%     c_M = basis( )

  end
  
end

