classdef ReNDSpace < ArraySpace
  %ReNDSpace
  %   Basic linear space of multidimensional arrays,
  %   whose elements may be represented by an array.
  %   
  % See also LinMatSpace
  
  methods 
    function this = ReNDSpace( dims )
      if numel(dims)==1
        % Detect if vector
        dims = [dims, 1];
      end
      this = this@ArraySpace( dims );
      
      % Build basis
      for k=1:this.vdim
        v_k = canvec( this.vdim, k );
        this.basis{k} = this.mat(v_k);
      end
    end
  end
  
end