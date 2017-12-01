classdef Re1DSpace < ArraySpace
  %Re1DSpace
  %   Most common linear space of Re^d vectors,
  %   whose elements may be represented by column vectors.
  % 
  % This class is rather an enriched wrapper around common
  % `vec` and `mat = vec^{-1}` operators.
  %   
  % See also LinMatSpace
  
  methods 
    function this = Re1DSpace( d )
      assert(numel(d)==1)
      dims = [d, 1]; % Column vectors
      this = this@ArraySpace( dims );
      
      % Build basis
      this.basis = num2cell( eye(this.dimVec), 1 );
    end
  end
  
end