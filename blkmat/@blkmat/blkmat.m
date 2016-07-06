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
% - The usual functionalities for matrix objects are extended.
%   For avoiding confusion, some of the traditional methods are
%   prepended the *blk* word, e.g. blknumel().

% Inspired by Kevin Murphy (www.cs.berkeley.edu/~murphyk), 21 October 1999
% This version by Jesus Briales, 3 May 2016

classdef blkmat
  
  properties
    rsizes, csizes
    rdict,  cdict
    storage
  end
  properties
    % Flags, set once at the constructor and never modified again
    is_rowRegular, is_colRegular
    is_labeled
  end
  
  methods   
    function this = blkmat( varargin )
      % (nrows, ncols, rsize, csize, M)
      % blkmat
      % Possible inputs are:
      % - NROWS, NCOLS, RSIZE, CSIZE
      % - [], [], RSIZES, CSIZES
      
      this.rdict = struct();
      this.cdict = struct();
      if nargin == 0
        % No arguments, default constructor
        this.rsizes = [];
        this.csizes = [];
        this.is_rowRegular = 0;
        this.is_colRegular = 0;
        this.storage = [];
        
      elseif nargin == 1
        % One single input argument
        if isa(varargin{1}, 'blkmat')
          % Copy constructor
          this = varargin{1};
          return % When clone copy is done, we can skip subsequent steps
        else
          % Build single block matrix from numeric input
          M = varargin{1};
          assert(isnumeric(M),'The input should be a numeric matrix');
          s = size(M);
          this.rsizes = s(1);
          this.csizes = s(2);
          this.storage = M;
        end
        
      else
        % There are more than 1 input argument
        if isa(varargin{1}, 'blkmat')
          % If 1st argument is blkmat, copy the input structure
          this = varargin{1};
          if nargin == 2
            % There is extra argument giving content initialization
            M = varargin{2};
          else
            M = 0; % Initialize to zero
          end
          
        else
          if nargin <= 3
            % The symmetric case, same dims and blk-sizes for rows and columns
            dims = varargin{1}; bsizes = varargin{2};
            [this.rsizes,this.rdict] = setupStructure( dims,bsizes );
            [this.csizes,this.cdict] = setupStructure( dims,bsizes );
          else
            % Set row structure
            rdims = varargin{1}; rsizes = varargin{3};
            [this.rsizes,this.rdict] = setupStructure( rdims,rsizes );
            % Set col structure
            cdims = varargin{2}; csizes = varargin{4};
            [this.csizes,this.cdict] = setupStructure( cdims,csizes );
          end
          
          % The extra argument for initialization can be the 3rd or 5th,
          % so it is the last if there is an odd number of arguments
          if mod(nargin,2)
            % If odd, there is extra argument giving content initialization
            M = varargin{end};
          else
            % If even, no initialization, set to 0
            M = 0;
          end
        end
        
        % Initialize storage field
        assert(isscalar(M) || all(size(M)==size(this)),...
          'blkmat: Wrong dim of the initialization matrix')
        if isscalar(M)
          this.storage = M*ones(size(this));
        else
          this.storage = M;
        end
      end
      
      % Set regularity flags
      this.is_rowRegular = (numel(unique(this.rsizes))==1);
      this.is_colRegular = (numel(unique(this.csizes))==1);
      % Set flag for labelled blkmat
      this.is_labeled = ~isempty(fieldnames(this.rdict)) || ...
                        ~isempty(fieldnames(this.cdict));
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
        if this.is_labeled
          S = map_lab2idx( this, S );
        end
        switch S.type
          case '()',
            [rows, cols] = extract_indices(this, S.subs);
            
            % If we reference a non-existent cell, expand the array if regular
            r = max(rows);
            if r > nrows(this)
              if this.is_rowRegular
                this.rsizes = repmat(this.rsizes(1), 1, r);
              else
                error('can''t expand row irregular blockmatrix');
              end
            end
            c = max(cols);
            if c > ncols(this)
              if this.is_colRegular
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
        
        if this.is_labeled % TODO: Change by is_labeled flag
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
        
    function [rows, cols] = extract_indices(this, subs)
      
      if length(subs)==1
        if ncols(this)==1 % col vector
          rows = subs{1};
          cols = 1;
        elseif nrows(this) == 1 % row vector
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
      if ischar(rows) && rows == ':'
        rows = 1:nrows(this);
      end
      
      if ischar(cols) && cols == ':'
        cols = 1:ncols(this);
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
    
    % Unary operators
    function B = uminus(A)
      B = A;
      B.storage = -A.storage;
    end
    function At = ctranspose(A)
      At = blkmat([],[],colsizes(A),rowsizes(A),A.storage');
    end
    function t = trace(A)
      t = trace(A.storage);
    end
    
    % Binary operators
    function C = plus(A,B)
      
      % Force blkmat matrices (in case one of the inputs is numeric)
      A = blkmat(A); B = blkmat(B);
      
      % Check consistent block dimensions
      assert( all(rowsizes(A)==rowsizes(B)) && ...
              all(colsizes(A)==colsizes(B)) )
      
      % Sum of internal storages
      temp = plain(A) + plain(B);
      
      % Store result in blkmat with the same structure
      C = blkmat([],[],rowsizes(A),colsizes(A),temp);
    end
    
    function C = minus(A,B)
      
      % Force blkmat matrices (in case one of the inputs is numeric)
      A = blkmat(A); B = blkmat(B);
      
      % Check consistent block dimensions
      assert( all(rowsizes(A)==rowsizes(B)) && ...
              all(colsizes(A)==colsizes(B)) )
      
      % Difference of internal storages
      temp = plain(A) + plain(B);
      
      % Store result in blkmat with the same structure
      C = blkmat([],[],rowsizes(A),colsizes(A),temp);
    end  
    
    function C = mtimes(A,B)
      % C = mtimes(A,B)
      % Some semantics in the product with block matrix:
      % - Product with a scalar simply scales all the content
      % - If the result is a single block (scalar or matrix)
      %   return raw Matlab matrix
      
      % Force blkmat matrices (in case one of the inputs is numeric)
      A = blkmat(A); B = blkmat(B);
      
      % Extract plain data
      matA = plain(A);
      matB = plain(B);
      
      % Valid cases are:
      % - A and B blkmat, with compatible size and blksize
      % - A or B scalar
      assert( numel(matA)==1 || numel(matB)==1 || ...
              all(colsizes(A)==rowsizes(B)) )
      
      temp = matA*matB;
      
      % Store result in blkmat with the resulting structure
      C = blkmat([],[],rowsizes(A),colsizes(B),temp);
      if numel(C) == 1
        % Result is a single block, return as simpler numeric matrix
        C = plain(C);
      end
    end
    
    function C = mldivide(A,B)
      % C = mldivide(A,B)
      
      % Force blkmat matrices (in case one of the inputs is numeric)
      A = blkmat(A); B = blkmat(B);
      
      % Extract plain data
      matA = plain(A);
      matB = plain(B);
      
      % Valid cases are:
      % - A and B blkmat, with compatible size and blksize
      % - A or B scalar
      assert( numel(matA)==1 || numel(matB)==1 || ...
              all(colsizes(A)==rowsizes(B)) )
      
      temp = matA\matB;
      
      % Store result in blkmat with the resulting structure
      C = blkmat([],[],rowsizes(A),colsizes(B),temp);
      if numel(C) == 1
        % Result is a single block, return as simpler numeric matrix
        C = plain(C);
      end
    end
    
    function C = mrdivide(A,B)
      % C = mldivide(A,B)
      
      % Force blkmat matrices (in case one of the inputs is numeric)
      A = blkmat(A); B = blkmat(B);
      
      % Extract plain data
      matA = plain(A);
      matB = plain(B);
      
      % Valid cases are:
      % - A and B blkmat, with compatible size and blksize
      % - A or B scalar
      assert( numel(matA)==1 || numel(matB)==1 || ...
              all(colsizes(A)==rowsizes(B)) )
      
      temp = matA/matB;
      
      % Store result in blkmat with the resulting structure
      C = blkmat([],[],rowsizes(A),colsizes(B),temp);
      if numel(C) == 1
        % Result is a single block, return as simpler numeric matrix
        C = plain(C);
      end
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
    function is = isregular(A), is = A.is_rowRegular && A.is_colRegular; end
    function l = labels(this)
      l = {cell2mat(fieldnames(this.rdict)),...
           cell2mat(fieldnames(this.cdict))};
    end
    
    function M = plain(this), M = this.storage; end

    function is = isempty(A)
      is = isempty(A.storage);
    end
    function E = blkisempty(A)
      % E = blkisempty(A)
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
    
    function disp(this)
      % Display the metadata about blk-matrix structure only
      % To display the matrix content, use plain method
      
      c = class(this.storage);
      if this.is_rowRegular && this.is_colRegular
        s = sprintf('%dx%d %s-matrix of %dx%d blocks',...
          nrows(this), ncols(this), c, rowsize(this), colsize(this));
      elseif this.is_rowRegular && ~this.is_colRegular
        s = sprintf('%dx- %s-matrix of %d x [%s] blocks',...
          nrows(this), c, rowsize(this), num2str(colsizes(this)));
      elseif ~this.is_rowRegular && this.is_colRegular
        s = sprintf('-x%d %s-matrix of [%s] x %d blocks',...
          ncols(this), c, num2str(rowsizes(this)), colsize(this));
      else
        s = sprintf('-x- %s-matrix of [%s] x [%s] blocks',...
          c, num2str(rowsizes(this)), num2str(colsizes(this)));
      end
      
      % Add information about labels, if the blkmat is labeled
      strLabels = [];
      if this.is_labeled
        strLabels = sprintf(', labeled as [%s] x [%s]',...
          num2str(cell2mat(fieldnames(this.rdict))),...
          num2str(cell2mat(fieldnames(this.cdict))) );
      end
      
      % Print to screen
      fprintf('\t%s%s\n',s,strLabels);
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