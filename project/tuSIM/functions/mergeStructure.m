% from https://nl.mathworks.com/matlabcentral/fileexchange/48231-structure-cat-merge-and-override

function A = mergeStructure(varargin)
%% MERGE_STRUCT Merge fields in scalar structures
%
% Syntax:
%     A = MERGE_STRUCT(A, B)
%     A = MERGE_STRUCT(A1, A2, A3, A4, ...)
%
% Input:
%     A,B...   - structures to merge
%
% Output:
%     A        - merged structure
%
% Comments:
%     Merges multiple scalar structures resulting in a single structure with the
%     fieldnames and values from all the input structures. 
%
% See also struct, fieldnames
%   Created by: Johan Winges
%   $Revision: 1.0$  $Date: 2014-10-21 16:00:00$
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