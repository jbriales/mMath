classdef Re1DSpace < LinMatSpace
  %Re1DSpace
  %   Most common linear space of Re^d vectors,
  %   whose elements may be represented by column vectors.
  % 
  % This class is rather an enriched wrapper around common
  % `vec` and `mat = vec^{-1}` operators.
  %   
  % See also LinMatSpace
  
  properties( Dependent )
    basis
  end
  
  methods 
    function this = Re1DSpace( d )
      assert(numel(d)==1)
      dims = [d, 1];
      this = this@LinMatSpace( dims );
      
      % Set dimension of the vector space
      this.vdim = d;
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
      % The canonical vectors are just the columns in the identity matrix
      c_E = num2cell( eye(this.dimVec), 1 );     
    end
  end
  
end