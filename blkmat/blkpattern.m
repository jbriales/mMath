classdef blkpattern
  %blkpattern Class that store the metadata in a block-matrix
  %   This class handles all the metadata in relation to a block-matrix.
  
  properties
    % Metadata
    sizes
    dict
    % Flags, set once at the constructor and never modified again
    is_regular
    is_labeled
  end
  properties (Dependent)
    labels
  end
  
  methods
    function this = blkpattern(varargin)
      if nargin == 0
        this.sizes = [];
        this.dict = struct();
        
      elseif nargin == 1
        if isa(varargin{1}, 'blkpattern')
          % Copy constructor
          this = varargin{1};
          return % When clone copy is done, we can skip subsequent steps
        end
        
      else
        % There are more than 1 input argument
        dims = varargin{1}; bsizes = varargin{2};
        [this.sizes,this.dict] = setupStructure( dims,bsizes );
      end
      
      % Set regularity flag
      this.is_regular = (numel(unique(this.sizes))==1);
      % Set flag for labelled pattern
      this.is_labeled = ~isempty(fieldnames(this.dict));
    end
    
    function l = get.labels(this)
      l = cell2mat(fieldnames(this.dict));
    end
    
    function is = eq(a,b)
      % Compare the two different blk patterns. If the patterns are
      % labeled, the labeling must be consistent for equality.
      is = all(a.sizes == b.sizes) && all(a.labels == b.labels);
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