classdef SkewSpace < LinMatSpace
  %SkewSpace
  %   Linear space of dxd skew-symmetric matrices.
  %   
  % See also LinMatSpace, SymSpace.
  
  methods 
    function this = SkewSpace( d )
      assert(numel(d)==1)
      this = this@LinMatSpace( [d,d] );
      
      % Compute corresponding dimension of vector space
      this.vdim = nchoosek(d,2);
      
      % Populate basis and index-correpondence tables
      % for Sym(d) vector space
      d = unique(this.dimMat);
      k=0;
      for i=1:d % traverse row-wise (this is symmetric!)
        for j=i+1:d
          k = k+1;
          % Get basis matrix element
          E_ij = skewPart(sparse(i,j,1,d,d,2));
%           c_E{k} = E_ij;
          this.basis{k} = full( E_ij );
          % Store subscripts
          this.m_ij(k,:) = [i,j];
          % Store linear index
          this.m_k(i,j) = k;
          this.m_k(j,i) = k;
        end
      end
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
    
  end
  
end