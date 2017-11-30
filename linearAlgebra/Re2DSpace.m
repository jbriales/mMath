classdef Re2DSpace < LinMatSpace
  %Re2DSpace([d1 d2]) or Re2DSpace(d) == Re2DSpace([d d])
  %   Linear space of Re^{m x d} rectangular matrices,
  %   whose elements may be represented by general matrices.
  %   
  % See also LinMatSpace, Re1DSpace.
  
  properties( Dependent )
    basis
  end
  
  methods 
    function this = Re2DSpace( dims )
      if numel(dims)==1
        % Detect if square
        dims = [dims, dims];
      end
      this = this@LinMatSpace( dims );
      
      this.vdim = prod(dims);
    end
    
    function v = vec(this, M)
      % Wrap standalone function
      v = vec(M);     
    end
    function M = mat(this, v)
      % Wrap standalone function
      M = mat(v,this.Mdim);
    end
    
    function c_E = get.basis( this )
      % Generate complete list of basis vectors for Re(d1xd2) vector space
      % The canonical vectors/matrices are every sparse matrix with (i,j)=1
      nrows = this.dimMat(1);
      ncols = this.dimMat(2);
      c_E = cell(1,this.dimVec);
      k=0;
      % NOTE: we traverse column-wise
      for j=1:ncols
        for i=1:nrows
          k = k+1;
          E_ij = sparse(i,j,1,nrows,ncols,1);
%           c_E{k} = E_ij;
          c_E{k} = full( E_ij );
        end
      end
    end
  end
  
end