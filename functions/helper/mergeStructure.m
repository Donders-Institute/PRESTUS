% from https://nl.mathworks.com/matlabcentral/fileexchange/48231-structure-cat-merge-and-override

function A = mergeStructure(varargin)
% MERGESTRUCTURE  Merge fields from multiple scalar structs into one
%
% Later structs overwrite fields with the same name from earlier structs.
% All inputs must be scalar structs. Adapted from Johan Winges (2014),
% MATLAB File Exchange 48231.
%
% Use as:
%   A = mergeStructure(A, B)
%   A = mergeStructure(A1, A2, A3, ...)
%
% Input:
%   varargin - two or more scalar structs to merge
%
% Output:
%   A - merged struct containing all fields from all inputs
%
% See also: STRUCT, FIELDNAMES, MERGESTRUCT

%% Merge structures:
% Find fieldnames of structures:
fAn   = cellfun(@(An) fieldnames(An), varargin,'un',0);
% Merge all structures with the first structure:
A = varargin{1};
if ~isscalar(A)
   error('merge_struct only works with scalar input structure A')
end
% Loop over all structures after the first:
for iA = 2:length(fAn)   
   % Loop over all fieldnames in each structure:
   for ifn = 1:length(fAn{iA})
      % Add the field from the structure to the first structure:
      if isscalar( varargin{iA} )
         % Structure elements are scalar:
         A.(fAn{iA}{ifn}) = varargin{iA}.(fAn{iA}{ifn});
      else
         % Structure to merge is an array, we position all data into a cell:
         tmpSize = size( varargin{iA} );
         tmpData = cell(tmpSize);
         [tmpData{:}] =  varargin{iA}.(fAn{iA}{ifn});
         A.(fAn{iA}{ifn}) = tmpData;         
      end
      % Note, this will override any fieldname in A if they exists in B.
   end
end