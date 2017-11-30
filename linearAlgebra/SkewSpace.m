classdef SkewSpace < LinMatSpace
  %SkewSpace
  %   Linear space of dxd skew-symmetric matrices.
  %   
  % See also LinMatSpace, SymSpace.
  
  properties( Dependent )
    basis
  end
  
  methods 
    function this = SkewSpace( d )
      assert(numel(d)==1)
      this = this@LinMatSpace( [d,d] );
      
      % Compute corresponding dimension of vector space
      this.vdim = nchoosek(d-1,2);
    end
    
    function v = vec(this, M)     
      error('TODO: Complete similarly to SymSpace')
    end
    function M = mat(this, v)
      error('TODO: Complete similarly to SymSpace')
    end
    function v = vecTr(M)
      error('TODO: Complete similarly to SymSpace')
    end
    function M = matTr(v)
      error('TODO: Complete similarly to SymSpace')
    end
    
    function c_E = get.basis( this )
      % Generate complete list of basis vectors for Skew(d) vector space
      d = unique(this.dimMat);
      c_E = cell(1,this.dimVec);
      k=0;
      for i=1:d
        for j=i+1:d
          k = k+1;
          E_ij = skewPart(sparse(i,j,1,d,d,2));
%           c_E{k} = E_ij;
          c_E{k} = full( E_ij );
        end
      end
    end
  end
  
end