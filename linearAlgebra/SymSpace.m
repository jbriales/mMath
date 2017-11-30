classdef SymSpace < LinMatSpace
  %SymSpace
  %   Linear space of dxd symmetric matrices.
  %   
  % See also LinMatSpace, ReSpace.
  
  properties( Dependent )
    basis
  end
  
  methods 
    function this = SymSpace( d )
      assert(numel(d)==1)
      this = this@LinMatSpace( [d,d] );
      
      % Compute corresponding dimension of vector space
      this.vdim = nchoosek(d,2);
    end
    
    function v = vec(this, M)     
      % Half vectorization for symmetric matrices
      % Take only elements in the lower triangular block of the matrix    
      bTrue = tril(true(this.dimMat));
      v = M( bTrue );
    end
    function M = mat(this, v)
      
      % Inverse of symmetric (half) vectorization for symmetric matrices
%       k = numel(x);
%       s = -0.5 + sqrt(0.25+2*k);
      s = this.dimMat;
      M = zeros(s)*(1); % Force same type as x
      
      % Insert vector to lower triangular block
      M( tril(true(s)) ) = v;
      
      % Scale the off-diagonal lower triangular block with 1/sqrt(2)
      M( ~triu(true(s)) ) = M( ~triu(true(s)) ) * 2;
      
      % Symmetrize matrix (x2 factor above to sum both terms)
      M = symmetrize( M );
    end
    
    function v = vecTr(M)
      % Symmetric (half) vectorization for symmetric matrices
      % This vectorization is scaled so that tr(ST)=vecs(S)'*vecs(T)
      % Taken from "Schäcke, K. (2013). On the kronecker product."
      
      s = this.dimMat;
      
      % Scale the off-diagonal lower triangular block with sqrt(2)
      M( ~triu(true(s)) ) = sqrt(2) * M( ~triu(true(s)) );
      
      % Extract lower triangular block to vector
      v = M( tril(true(s)) );
    end
    
    function M = matTr(v)
      % Inverse of symmetric (half) vectorization for symmetric matrices
      % Scaled half vectorization so that tr(ST)=vecs(S)'*vecs(T)
      % Taken from "Schäcke, K. (2013). On the kronecker product."
      
%       k = numel(v);
%       s = -0.5 + sqrt(0.25+2*k);
      s = this.dimMat;
      M = zeros(s);
      
      % Insert vector to lower triangular block
      M( tril(true(s)) ) = v;
      
      % Scale the off-diagonal lower triangular block with 1/sqrt(2)
      M( ~triu(true(s)) ) = M( ~triu(true(s)) ) * 2 / sqrt(2);
      
      % Symmetrize matrix (x2 factor above to sum both terms
      M = symmetrize( M );
      
    end
    
    function c_E = get.basis( this )
      % Generate complete list of basis vectors for Sym(d) vector space
      d = unique(this.dimMat);
      c_E = cell(1,this.dimVec);
      k=0;
      for i=1:d
        for j=i:d
          k = k+1;
          E_ij = symPart(sparse(i,j,1,d,d,2));
%           c_E{k} = E_ij;
          c_E{k} = full( E_ij );
        end
      end
    end
  end
  
end