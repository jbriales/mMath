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
%
% Notes:
% - Once the blkmat object has been built, its structure is fixed.
%   Only the content can be modified via subscripted assignment.

% Inspired by Kevin Murphy (www.cs.berkeley.edu/~murphyk), 21 October 1999
% This version by Jesus Briales, 3 May 2016

classdef blkmat
  %BLKMAT Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    rsizes, csizes
    rdict, cdict
    storage
  end
  properties
    % Flags, set once at the constructor and never modified again
    row_regular, col_regular
  end
  properties
    isLabeled
  end
  
  methods   
    function this = blkmat( varargin )
      % (nrows, ncols, rsize, csize, M)
      % blkmat
      % Possible inputs are:
      % - NROWS, NCOLS, RSIZE, CSIZE
      % - [], [], RSIZES, CSIZES
%       keyboard
      
      this.rdict = struct();
      this.cdict = struct();
      if nargin == 0 % default constructor
        this.rsizes = [];
        this.csizes = [];
        this.row_regular = 0;
        this.col_regular = 0;
        this.storage = [];
        
      elseif isa(varargin{1}, 'blkmat')
        % If 1st argument is blkmat, this is a copy constructor
        this = varargin{1};
        if nargin == 2
          % There is extra argument giving content initialization
          M = varargin{2};
        else
          M = 0; % Initialize to zero
        end
        
      else
        % Set row structure
        [this.rsizes,this.rdict] = setupStructure( varargin{[1 3]} );
        this.row_regular = (numel(unique(this.rsizes))==1);
        
        % Set col structure
        [this.csizes,this.cdict] = setupStructure( varargin{[2 4]} );
        this.col_regular = (numel(unique(this.csizes))==1);
        
        if nargin == 5
          % There is extra argument giving content initialization
          M = varargin{5};
        else
          M = 0;
        end
      end
      
      % Set flag for labelled blkmat
      this.isLabeled = ~isempty(fieldnames(this.rdict)) || ...
                        ~isempty(fieldnames(this.cdict));
      
      % Initialize storage field
      assert(isscalar(M) || all(size(M)==size(this)),...
        'blkmat: Wrong dim of the initialization matrix')
      if isscalar(M)
        this.storage = M*ones(size(this));
      else
        this.storage = M;
      end
    end
       
    % Set the number of arguments in subscript operations
    % Avoid conflicts with (overloaded) numel method
    function n = numArgumentsFromSubscript(this,~,callingContext)
      switch callingContext
        case matlab.mixin.util.IndexingContext.Statement
          n = 1; % nargout for indexed reference used as statement
        case matlab.mixin.util.IndexingContext.Expression
          n = 1; % nargout for indexed reference used as function argument
        case matlab.mixin.util.IndexingContext.Assignment
          n = 1; % nargin for indexed assignment
      end
    end
    
    function S = map_lab2idx( this, S )
      if S.type == '.'
        % If field accessor, convert into usual indexing
        % (only one-char index per dimension can be used)
        S.type = '()';
        S.subs = num2cell(S.subs);
      end
      if strcmp(S.type,'()')
        subs = S.subs;
        switch length(subs)
          case 2
            % Row and column indeces are given
            S.subs{1} = label2idx( this.rdict, S.subs{1} );
            S.subs{2} = label2idx( this.cdict, S.subs{2} );
          case 1
            % Case a single index is given (blk-vector)
            if ncols(this)==1 % col vector
              S.subs{1} = label2idx( this.rdict, S.subs{1} );
            elseif nrows(this) == 1 % row vector
              S.subs{1} = label2idx( this.cdict, S.subs{1} );
            else
              error('This is not a blk-vector, specify two indices for matrices');
            end
          otherwise
            error('Blk-mat has maximum 2 dimensions');
        end

      end
      
      function idxs = label2idx( dict, labels )
        if ~strcmp(labels,':')
          % Check all the labels are contained in the dictionary
          labelsNotInDict = cell2mat( setdiff(labels,fieldnames(dict)) );
          assert(isempty(labelsNotInDict),...
                 'Labels %s don''t exist in the matrix dimension',labelsNotInDict);
          % Read indices corresponding to labels from the dictionary
          n = numel(labels);
          idxs = zeros(1,n);
          for j=1:n
            l = labels(j);
            idxs(j) = dict.(l);
          end
        else
          idxs = ':';
        end
      end
    end
    
    function this = subsasgn(this,S,B)
      
      if length(S)==1 % A(i,j)
        if this.isLabeled
          S = map_lab2idx( this, S );
        end
        switch S.type
          case '()',
            [rows, cols] = extract_indices(this, S.subs);
            
            % If we reference a non-existent cell, expand the array if regular
            r = max(rows);
            if r > nrows(this)
              if this.row_regular
                this.rsizes = repmat(this.rsizes(1), 1, r);
              else
                error('can''t expand row irregular blockmatrix');
              end
            end
            c = max(cols);
            if c > ncols(this)
              if this.col_regular
                this.csizes = repmat(this.csizes(1), 1, c);
              else
                error('can''t expand column irregular blockmatrix');
              end
            end
            
            scalar_rows = blk2sub(rows, this.rsizes);
            scalar_cols = blk2sub(cols, this.csizes);
            this.storage(scalar_rows, scalar_cols) = B;
          otherwise,
            error(['Unrecognized subscript operator ' S.type]);
        end
      else
        error(['Unrecognized subscript operator ' S]);
      end
    end
    
    function varargout = subsref(this,S)
            
      if length(S)==1 % A(i,j)
        
        if this.isLabeled % TODO: Change by isLabeled flag
          S = map_lab2idx( this, S );
        end
        
        switch S.type
          case '()',
            % Use numerical indeces for blocks in the matrix
            [rows, cols] = extract_indices(this, S.subs);
            if max(rows) > nrows(this) | max(cols) > ncols(this)
              error(['index ' num2str(rows) ', ' num2str(cols) ' is out of bounds']);
            end
            varargout{1} = this.storage(blk2sub(rows, rowsizes(this)), blk2sub(cols, colsizes(this)));
            
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
    function l = labels(this)
      l = {cell2mat(fieldnames(this.rdict)),...
           cell2mat(fieldnames(this.cdict))};
    end
    
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
    
    function disp(A)
      % To display the matrix content, use plain method
      if A.row_regular && A.col_regular
        s = sprintf('%dx%d %s-matrix of %dx%d blocks',...
          nrows(A), ncols(A), class(A.storage), rowsize(A), colsize(A));
      elseif A.row_regular && ~A.col_regular
        s = sprintf('%dx- %s-matrix of %d x [%s] blocks',...
          nrows(A), class(A.storage), rowsize(A), num2str(colsizes(A)));
      elseif ~A.row_regular && A.col_regular
        s = sprintf('-x%d %s-matrix of [%s] x %d blocks',...
          ncols(A), class(A.storage), num2str(rowsizes(A)), colsize(A));
      else
        s = sprintf('-x- %s-matrix of [%s] x [%s] blocks',...
          class(A.storage), num2str(rowsizes(A)), num2str(colsizes(A)));
      end
      fprintf('\t%s\n',s);
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

function [sizes,dict] = setupStructure( arg1, sizes )
% Preallocate dictionary
dict = struct();

if ischar(arg1)
  rtags = arg1;
  nblks = numel(rtags);
  for i = 1:nblks
    l = rtags(i);
    dict.(l) = i;
  end
else
  nblks = arg1;
  if isempty(nblks)
    % If empty nblks, this is take from dimension of rsizes
    nblks = numel(sizes);
  end
end
if isscalar(sizes)
  % If we are given a single block-size, repeat this
  sizes = repmat(sizes, 1, nblks);
end
assert(nblks==numel(sizes),...
  'Number of blocks (%d) and list of block dims (%d) should be consistent',...
  nblks, numel(sizes));
end