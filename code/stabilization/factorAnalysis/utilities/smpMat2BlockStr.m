function blockStr = smpMat2BlockStr(x)
% This is a function to convert a matrix of samples with missing values
% indicated by NAN to a structure organized by block, where each block has
% the same non-missing variables.
%
% Usage: blockStr = smpMat2BlockStr(x)
%
% Inputs: 
%
%   x - A matrix of samples.  Rows are samples, columns are variables.
%   Missing values are indicated by NAN. 
%
% Outputs: 
%
%   blockStr - A structure of length B, where B is the number of blocks.
%   Each entry will have the follwing fields:
%
%       nSmps - The number of samples that are in this block. 
%
%       blockXInds - A vector listing the columns in the x matrix for the
%       variables that are represented in this block - that is the columns
%       that were not NAN values in the original x matrix. 
%
%       origRows - a vector listing the rows the samples of x the samples
%       for this block were taken from. 
%
%       obsX - A matrix of size N by D, where N is the number of samples in
%       the block and D is the number of dimensions that are non NAN for
%       this block (will correspond to the length of blockXInds).  The
%       formula obsX = x(blockStr(b).origRows, blockStr(b).blockXInds) 
%       relates the samples in obsX to the original data. 
%      
% Author: William Bishop, bishopw@janelia.hhmi.org

% See which values in x are NAN
xNotNanInds = ~isnan(x);

% Now we get unique blocks
[uniqueRows, ~, rowInds] = unique(xNotNanInds, 'rows');

nBlocks = size(uniqueRows,1);

for b = 1:nBlocks
    origRows = find(rowInds == b);
    blockXInds = find(uniqueRows(b,:));
    
    blockStr(b).nSmps = length(origRows);
    blockStr(b).blockXInds = blockXInds; 
    blockStr(b).origRows = origRows; 
    blockStr(b).obsX = x(origRows, blockXInds); 
end