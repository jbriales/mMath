classdef ReNDSpace < LinMatSpace
  %ReNDSpace
  %   Basic linear space of multidimensional arrays,
  %   whose elements may be represented by an array.
  %   
  % See also LinMatSpace
  
  properties( Dependent )
    basis
  end
  
  methods 
    function this = ReNDSpace( dims )
      if numel(dims)==1
        % Detect if vector
        dims = [dims, 1];
      end
      this = this@LinMatSpace( dims );
      
      this.vdim = prod(dims);
    end
    
    function v = vec(this, M)
      
      % Vectorize column-wise
      v = M(:);
      
    end
    function M = mat(this, v)
      
      % Reshape vector into original array dimensions
      M = reshape(v,this.Mdim);
      
    end
    
    function c_M = get.basis( this )
      
      c_M = cell(1,this.vdim);
      for k=1:this.vdim
        v_k = sparse(k,1,1,this.vdim,1,1);
        c_M{k} = this.mat(v_k);
      end
      
    end
  end
  
end