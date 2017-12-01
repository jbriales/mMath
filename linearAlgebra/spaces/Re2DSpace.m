classdef Re2DSpace < ArraySpace
  %Re2DSpace([d1 d2]) or Re2DSpace(d) == Re2DSpace([d d])
  %   Linear space of Re^{m x d} rectangular matrices,
  %   whose elements may be represented by general matrices.
  %   
  % See also LinMatSpace, Re1DSpace.
  
  methods 
    function this = Re2DSpace( dims )
      if numel(dims)==1
        % Detect if square
        dims = [dims, dims];
      end
      this = this@ArraySpace( dims );
      
      % Generate complete list of basis vectors for Re(d1xd2) vector space
      % The canonical vectors/matrices are every sparse matrix with (i,j)=1
      nrows = this.dimMat(1);
      ncols = this.dimMat(2);
      k=0;
      % NOTE: we traverse column-wise
      for j=1:ncols
        for i=1:nrows
          k = k+1;
          E_ij = sparse(i,j,1,nrows,ncols,1);
%           c_E{k} = E_ij;
          this.basis{k} = full( E_ij );
        end
      end
    end    
  end
  
end