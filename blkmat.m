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
    M
  end
  
  methods   
    function this = blkmat(nrows, ncols, rsize, csize, M)
      if nargin == 0 % default constructor
        this.rsizes = [];
        this.csizes = [];
        this.row_regular = 0;
        this.col_regular = 0;
        this.M = [];
        
      elseif isa(nrows, 'blkmat')
        this = nrows; % identity function
        
      else
        if isempty(nrows)
          this.rsizes = rsize;
          this.row_regular = 0;
        else
          this.rsizes = repmat(rsize, 1, nrows);
          this.row_regular = 1;
        end
        
        if isempty(ncols)
          this.csizes = csize;
          this.col_regular = 0;
        else
          this.csizes = repmat(csize, 1, ncols);
          this.col_regular = 1;
        end
        
        if nargin < 5
          this.M = zeros(sum(this.rsizes), sum(this.csizes));
        else
          % Check initialization is scalar or dimensions are coherent
          assert(isscalar(M) || all(size(M)==size(this)))
          if isscalar(M)
            this.M = M*ones(size(this));
          else
            this.M = M;
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
            A.M(scalar_rows, scalar_cols) = B;
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
            B = A.M(blk2sub(rows, rowsizes(A)), blk2sub(cols, colsizes(A)));
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
      % Check consistent block dimensions
%       assert(all(rowsizes(A)==rowsizes(B))&&all(colsizes(A)==colsizes(B)))
      % Get sum of entire content
%       temp = A.M + B.M;
      if isa(A,'blkmat'), M_A = A.M; else M_A = A; end
      if isa(B,'blkmat'), M_B = B.M; else M_B = B; end
      temp = M_A + M_B;
      
      % TODO: Create according to regular block-matrix or not
      obj = getobj(A,B);
      if ~obj.row_regular && ~obj.col_regular
        C = blkmat([],[],rowsizes(obj),colsizes(obj),temp);
      else
        error('TODO yet');
      end
    end
    
    function C = minus(A,B)
      % Check consistent block dimensions
%       assert(all(rowsizes(A)==rowsizes(B))&&all(colsizes(A)==colsizes(B)))
      % Get sum of entire content
%       temp = A.M + B.M;
      if isa(A,'blkmat'), M_A = A.M; else M_A = A; end
      if isa(B,'blkmat'), M_B = B.M; else M_B = B; end
      temp = M_A - M_B;
      
      % TODO: Create according to regular block-matrix or not
      obj = getobj(A,B);
      if ~obj.row_regular && ~obj.col_regular
        C = blkmat([],[],rowsizes(obj),colsizes(obj),temp);
      else
        error('TODO yet');
      end
    end
    
    function C = mtimes(A,B)
      % Check consistent block dimensions
      % TODO
%       assert(all(rowsizes(A)==rowsizes(B))&&all(colsizes(A)==colsizes(B)))
      % Get sum of entire content
      if isa(A,'blkmat'), M_A = A.M; obj=A; else M_A = A; end
      if isa(B,'blkmat'), M_B = B.M; obj=B; else M_B = B; end
      temp = M_A * M_B;
      
      % TODO: Create according to regular block-matrix or not
      if ~obj.row_regular && ~obj.col_regular
        C = blkmat([],[],rowsizes(obj),colsizes(obj),temp);
      else
        error('TODO yet');
      end
    end
    
    function At = ctranspose(A)
      % Transpose internal values of temporary A
      A.M = A.M';
      % Create new object
      At = blkmat(A);
    end
    
    function r = rowsize(A),  r = A.rsizes(1); end
    function r = rowsizes(A), r = A.rsizes; end
    function r = nrows(A),    r = length(A.rsizes); end
    function c = colsize(A),  c = A.csizes(1); end
    function c = colsizes(A), c = A.csizes; end
    function c = ncols(A),    c = length(A.csizes); end
    function s = size(A)
      s = [sum(A.rsizes), sum(A.csizes)];
    end
    
    function E = isempty(A)
    % E = isempty(A)
    % Check for each block if all elements are zero
    nr = nrows(A); nc = ncols(A);
    E = false(nr,nc);
    for i=1:nr
      for j=1:nc
        B = A.M(blk2sub(i,rowsizes(A)), blk2sub(j,colsizes(A)));
        E(i,j) = all(all(B==0));
      end
    end
    end
    
    function meta(A)
    % meta(A)
    % Display metadata for the block matrix
      if A.row_regular & A.col_regular
        s = sprintf('%dx%d %s-matrix of %dx%d blocks',...
          nrows(A), ncols(A), class(A.M), rowsize(A), colsize(A));
      elseif A.row_regular & ~A.col_regular
        s = sprintf('%dx- %s-matrix of %d x [%s] blocks',...
          nrows(A), class(A.M), rowsize(A), num2str(colsizes(A)));
      elseif ~A.row_regular & A.col_regular
        s = sprintf('-x%d %s-matrix of [%s] x %d blocks',...
          ncols(A), class(A.M), num2str(rowsizes(A)), colsize(A));
      else
        s = sprintf('-x- %s-matrix of [%s] x [%s] blocks',...
          class(A.M), num2str(rowsizes(A)), num2str(colsizes(A)));
      end
      disp(' ');
      disp([' ' inputname(1) ' = ' s ]);
      disp(' ');
    end
    
    function disp(A)
      meta(A)
      % TODO: Show only if dimensions is below some threshold
      disp(A.M)
    end
    
    function spy(A)
      spy(A.M);
      % If exists spyblck, add lines between blocks
      % TODO: Check function exists
      plt.spyblk(rowsizes(A),colsizes(A));
    end
    
  end

end