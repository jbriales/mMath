% BLKMAT Blockmatrix class.
% This class works as a convenient wrapper to improve code readability
% when working with matrices that are conceptually distributed blockwise.
%
% M = BLKMAT(NROWS, NCOLS, RSIZE, CSIZE) creates an NROWS by NCOLS blockmatrix,
% where the size of each block is RSIZE by CSIZE. This is called a regular blockmatrix.
%
% It is also possible for the blocks to have different sizes along a dimension.
% M = BLKMAT([], NCOLS, RSIZES, CSIZE) creates a NROWS by NCOLS blockmatrix,
% where NROWS = LENGTH(RSIZES), and the size of block (i,j) is RSIZES(I) by CSIZE.
% This is called a row irregular, column regular blockmatrix.
% M = BLKMAT(NROWS, [], RSIZE, CSIZES) creates a NROWS by NCOLS blockmatrix,
% where NCOLS = LENGTH(CSIZES), and the size of block (i,j) is RSIZE by CSIZES(J)
% This is called a row regular, column irregular blockmatrix.
% M = BLKMAT([], [], RSIZES, CSIZES) creates a NROWS by NCOLS blockmatrix,
% where the size of block (i,j) is RSIZES(I) by CSIZES(J).
% This is called a fully irregular matrix.
%
% Note that M = BLKMAT(NROWS, NCOLS, RSIZE, CSIZE) is equivalent to
% M = BLKMAT([], [], RSIZE*ONES(1,NROWS), CSIZE*ONES(1,NCOLS)),
% except that the latter is considered irregular.
%
% Use 'x = M(i,j)' to assign block (i,j) to x.
% If this block does not exist, an error will be raised.
% If ncols(M)=1, it suffices to write 'x = M(i)', and similarly for row vectors.
% You can also use ranges, e.g., M(:,:) or M(:) or M(1:3)
%
% Use 'M(i,j) = x' to assign x to block (i,j).
% If this block does not exist, it will be created,
% providing M is regular in the out-of-bounds dimension.
% (If M is irregular, an error will be raised.)
%
% The following functions are defined in the obvious way:
% NROWS(M), NCOLS(M)
% ROWSIZES(M), COLSIZES(M) if irregular
% ROWSIZE(M), COLSIZE(M) if regular
%
% No other operations can be performed on blockmatrices.
% To do arithmetic, you need to operate on the underlying scalar matrix,
% which can be accessed using M(:,:).
%
% EXAMPLES
%   M = blockmatrix(2, 2, 2, 2)
%   M(1,1) = eye(2); M(2,2) = -eye(2).
%
%     1     0     0     0
%     0     1     0     0
%     0     0    -1     0
%     0     0     0    -1
%
%   M = blockmatrix([], [], [2 3], [2 1])
%   M(:,:) = rand(5, 3) looks like this:
%
%    x x   x
%    x x   x
%
%    x x   x
%    x x   x
%    x x   x

% Written by Kevin Murphy (www.cs.berkeley.edu/~murphyk), 21 October 1999
% Refactored by Jesus Briales, 3 May 2016

classdef blkmat
  %BLKMAT Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    rsizes, csizes
    row_regular, col_regular
    storage
  end
  
  methods   
    function this = blkmat(nrows, ncols, rsize, csize, M)
      if nargin == 0 % default constructor
        this.rsizes = [];
        this.csizes = [];
        this.row_regular = 0;
        this.col_regular = 0;
        this.storage = [];
        
      elseif isa(nrows, 'blkmat')
        this = nrows; % identity function
        
      else
        if isempty(nrows)
          this.rsizes = rsize;
        else
          this.rsizes = repmat(rsize, 1, nrows);
        end
        if numel(unique(this.rsizes))==1
          this.row_regular = 1;
        else
          this.row_regular = 0;
        end
        
        if isempty(ncols)
          this.csizes = csize;
        else
          this.csizes = repmat(csize, 1, ncols);
        end
        if numel(unique(this.csizes))==1
          this.col_regular = 1;
        else
          this.col_regular = 0;
        end
        
        if nargin < 5
          this.storage = zeros(sum(this.rsizes), sum(this.csizes));
        else
          % Check initialization is scalar or dimensions are coherent
          assert(isscalar(M) || all(size(M)==size(this)),...
            'Wrong dimensions of the initialization matrix')
          if isscalar(M)
            this.storage = M*ones(size(this));
          else
            this.storage = M;
          end
        end
      end
    end
    
    function A = subsasgn(A,S,B)
      
      if length(S)==1 % A(i,j)
        switch S.type
          case '()',
            [rows, cols] = extract_indices(A, S.subs);
            
            % If we reference a non-existent cell, expand the array if regular
            r = max(rows);
            if r > nrows(A)
              if A.row_regular
                A.rsizes = repmat(A.rsizes(1), 1, r);
              else
                error('can''t expand row irregular blockmatrix');
              end
            end
            c = max(cols);
            if c > ncols(A)
              if A.col_regular
                A.csizes = repmat(A.csizes(1), 1, c);
              else
                error('can''t expand column irregular blockmatrix');
              end
            end
            
            scalar_rows = blk2sub(rows, A.rsizes);
            scalar_cols = blk2sub(cols, A.csizes);
            A.storage(scalar_rows, scalar_cols) = B;
          otherwise,
            error(['Unrecognized subscript operator ' S.type]);
        end
      else
        error(['Unrecognized subscript operator ' S]);
      end
    end
    
    function B = subsref(A,S)
      
      if length(S)==1 % A(i,j)
        switch S.type
          case '()',
            [rows, cols] = extract_indices(A, S.subs);
            if max(rows) > nrows(A) | max(cols) > ncols(A)
              error(['index ' num2str(rows) ', ' num2str(cols) ' is out of bounds']);
            end
            B = A.storage(blk2sub(rows, rowsizes(A)), blk2sub(cols, colsizes(A)));
          otherwise,
            error(['Unrecognized subscript operator ' S.type]);
        end
      else
        error(['Unrecognized subscript operator ' S]);
      end
    end
    
    function [rows, cols] = extract_indices(A, subs)
      
      if length(subs)==1
        if ncols(A)==1 % col vector
          rows = subs{1};
          cols = 1;
        elseif nrows(A) == 1 % row vector
          cols = subs{1};
          rows = 1;
        else
          error('must specify two indices');
        end
      else
        rows = subs{1};
        cols = subs{2};
      end
      
      % just using the test "rows == ':'" produces strange results if rows is 58, since char(58)==':'
      if isstr(rows) & rows == ':'
        rows = 1:nrows(A);
      end
      
      if isstr(cols) & cols == ':'
        cols = 1:ncols(A);
      end
    end
    
    function obj = getobj(varargin)
      % obj = getobj(varargin)
      % Returns the first found object of type blkmat
      for i=1:numel(varargin)
        if isa(varargin{i},'blkmat')
          obj = varargin{i};
          return
        end
      end
    end
    
    function C = plus(A,B)
      
      % If any input is not blkmat, convert to blkmat with single block
      if ~isa(A,'blkmat')
        A = blkmat(1,1,size(A,1),size(A,2),A);
      end
      if ~isa(B,'blkmat')
        B = blkmat(1,1,size(B,1),size(B,2),B);
      end
      
      % Check consistent block dimensions
      assert( all(rowsizes(A)==rowsizes(B)) && ...
              all(colsizes(A)==colsizes(B)) )
      
      % Extract plain data
      matA = plain(A); matB = plain(B);

      % Perform sum in internal storage
      temp = matA + matB;
      
      % Store result in blkmat with the same structure
      C = blkmat([],[],rowsizes(A),colsizes(A),temp);
    end
    
    function C = minus(A,B)
      
      % If any input is not blkmat, convert to blkmat with single block
      if ~isa(A,'blkmat')
        A = blkmat(1,1,size(A,1),size(A,2),A);
      end
      if ~isa(B,'blkmat')
        B = blkmat(1,1,size(B,1),size(B,2),B);
      end
      
      % Check consistent block dimensions
      assert( all(rowsizes(A)==rowsizes(B)) && ...
              all(colsizes(A)==colsizes(B)) )
      
      % Extract plain data
      matA = plain(A); matB = plain(B);

      % Perform sum in internal storage
      temp = matA - matB;
      
      % Store result in blkmat with the same structure
      C = blkmat([],[],rowsizes(A),colsizes(A),temp);
    end
    
    function B = uminus(A)
      B = A;
      B.storage = -A.storage;
    end      
    
    function C = mtimes(A,B)
      % C = mtimes(A,B)
      % Blk-matrix product
      % Some semantics in the product with block matrix:
      % - The block-matrix should be regular?
      % - The block-wise product must be consistent
      % - Product with a scalar simply scales all the content
      % - If the result is a single block (scalar or matrix)
      %   return raw Matlab matrix
      
      % If any input is not blkmat, convert to blkmat with single block
      if ~isa(A,'blkmat')
        A = blkmat(1,1,size(A,1),size(A,2),A);
      end
      if ~isa(B,'blkmat')
        B = blkmat(1,1,size(B,1),size(B,2),B);
      end
      
      % Extract plain data
      matA = plain(A);
      matB = plain(B);
      
      % Valid cases are:
      % - A and B blkmat, with compatible size and blksize
      % - A or B scalar
      assert( numel(matA)==1 || numel(matB)==1 || ...
              (ncols(A)==nrows(B) && colsize(A)==rowsize(B)) )
      
      temp = matA*matB;
      % Build the corresponding blkmat for the output
      if numel(matA)==1
        nr=nrows(B); nc=ncols(B); rs=rowsize(B); cs=colsize(B);
      elseif numel(matB)==1
        nr=nrows(A); nc=ncols(A); rs=rowsize(A); cs=colsize(A);
      else
        nr=nrows(A); nc=ncols(B); rs=rowsize(A); cs=colsize(B);
      end
      C = blkmat(nr,nc,rs,cs,temp);
      if numel(C) == 1
        % Scalar result, return simple number
        C = plain(C);
      end
    end
    
    function C = mldivide(A,B)
      % Check consistent block dimensions
      % TODO
%       assert(all(rowsizes(A)==rowsizes(B))&&all(colsizes(A)==colsizes(B)))
      % Get sum of entire content
      
      assert( isa(A,'blkmat') && isa(B,'blkmat') );
      assert( ncols(A)==nrows(B) )
      % TODO: Assert compatible block size
      if isa(A,'blkmat'), M_A = A.storage; obj=A; else M_A = A; end
      if isa(B,'blkmat'), M_B = B.storage; obj=B; else M_B = B; end
      temp = A.storage \ B.storage;
      
      % TODO: Create according to regular block-matrix or not
      if ~obj.row_regular && ~obj.col_regular
        C = blkmat([],[],rowsizes(obj),colsizes(obj),temp);
      else
        error('TODO yet');
      end
    end
    
    function At = ctranspose(A)
      % Create new transposed object
      At = blkmat(ncols(A),nrows(A),colsize(A),rowsize(A),A.storage');
      % TODO: Change if not regular
    end
    
    function t = trace(A)
      t = trace(A.storage);
    end
    
    function r = rowsize(A),  r = A.rsizes(1); end
    function r = rowsizes(A), r = A.rsizes; end
    function r = nrows(A),    r = length(A.rsizes); end
    function c = colsize(A),  c = A.csizes(1); end
    function c = colsizes(A), c = A.csizes; end
    function c = ncols(A),    c = length(A.csizes); end
    function s = size(A)
      % TODO: Refactor so that size counts nrows and ncols
      % For complete size use size(A(:,:)) instead
      s = [sum(A.rsizes), sum(A.csizes)];
    end
    function n = numel(A), n = nrows(A)*ncols(A); end
    function s = blksize(A)
      % NOTE: regular blkmat should be asserted externally
      s = [rowsize(A),colsize(A)];
    end
    function isit = isregular(A), isit = A.row_regular && A.col_regular; end
    
    function M = plain(this), M = this.storage; end

    function e = isempty(A)
      e = isempty(A.storage);
    end
    function E = blkisempty(A)
      % E = isempty(A)
      % Check for each block if all elements are zero
      nr = nrows(A); nc = ncols(A);
      E = false(nr,nc);
      for i=1:nr
        for j=1:nc
          B = A.storage(blk2sub(i,rowsizes(A)), blk2sub(j,colsizes(A)));
          E(i,j) = all(all(B==0));
        end
      end
    end
    
    function meta(A)
    % meta(A)
    % Display metadata for the block matrix
      if A.row_regular & A.col_regular
        s = sprintf('%dx%d %s-matrix of %dx%d blocks',...
          nrows(A), ncols(A), class(A.storage), rowsize(A), colsize(A));
      elseif A.row_regular & ~A.col_regular
        s = sprintf('%dx- %s-matrix of %d x [%s] blocks',...
          nrows(A), class(A.storage), rowsize(A), num2str(colsizes(A)));
      elseif ~A.row_regular & A.col_regular
        s = sprintf('-x%d %s-matrix of [%s] x %d blocks',...
          ncols(A), class(A.storage), num2str(rowsizes(A)), colsize(A));
      else
        s = sprintf('-x- %s-matrix of [%s] x [%s] blocks',...
          class(A.storage), num2str(rowsizes(A)), num2str(colsizes(A)));
      end
      disp(' ');
      disp([' ' inputname(1) ' = ' s ]);
      disp(' ');
    end
    
    function disp(A)
      meta(A)
      % TODO: Show only if dimensions is below some threshold
      disp(A.storage)
    end
    
    function spy(A)
      spy(A.storage);
      % If exists spyblck, add lines between blocks
      % TODO: Check function exists
      plt.spyblk(rowsizes(A),colsizes(A));
    end
    
    function out = full(A)
      out = full(A.storage);
    end
    
  end

end