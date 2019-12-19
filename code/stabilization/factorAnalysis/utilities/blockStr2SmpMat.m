function x = blockStr2SmpMat(blockStr)
% This is a function to convert data stored in a block structure to
% matrix with missing values indicated by NAN
%
% Usage: x = blockStr2SmpMat(blockStr)
%
% Inputs: 
%
%   blockStr - data stored in a block structure.  See the function
%   smpMat2BlockStr for the definition of this structure. 
%
% Outputs: 
%
%   x - data in matrix form.  Rows are samples and variables are columns.
%   Missing values are indicated by NAN. 
%
% Author: William Bishop, bishopw@janelia.hhmi.org

nBlocks = length(blockStr);

nBlockSmps = [blockStr.nSmps];
nTotalSmps = sum(nBlockSmps); 

nCols = max(cellfun(@(x) max(x), {blockStr.blockXInds}));

x = nan(nTotalSmps, nCols); 
for b = 1:nBlocks
    x(blockStr(b).origRows, blockStr(b).blockXInds) = blockStr(b).obsX; 
end