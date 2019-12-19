function alignRows = identifyStableLoadingRows(m1, m2, nStableRows, th)
% Identifies stable rows of two loading matrices, which are not expressed
% in the same basis.
%
% This algorithm works by iteratively tring to align the two loading
% matrices and then identifying and removing rows which are the most
% different after each alignment.
%
% Usage: alignRows = identifyStableLoadingRows(m1, m2, nStableRows, th, varargin)
%
% Inputs:
%
%   m1, m2: Two loading matrices to compare
%
%   nStableRows: the number of rows to use for alignment
%
%   th: the threshold to use when screening out rows for possibly
%   alignment.  Any row which has an l_2 norm less than th in either m_1 or
%   m_2 will not be considered when selecting rows to use for alignment.
%
% Outputs:
%
%   alignRows: A boolean array of length equal to the number of rows in the
%   loading matrices, with values of true indicating corresponding stable
%   rows in m1 and m2. 
%
% Author: William Bishop, bishopw@janelia.hhmi.org

% Identify rows of m1 and m2 that are too small
m1Norms = sqrt(sum(m1.^2,2));
smallM1Rows = m1Norms < th;

m2Norms = sqrt(sum(m2.^2,2));
smallM2Rows = m2Norms < th;

smallRows = smallM1Rows | smallM2Rows;
cleanRows = find(~smallRows);
nCleanRows = length(cleanRows);

m1Clean = m1(cleanRows,:);
m2Clean = m2(cleanRows,:);


curRows = 1:nCleanRows;
nDropRows = max(nCleanRows - nStableRows, 0);
for i = 1:nDropRows
    % Perform alignment using the rows that were identified last iteration
    curM1 = m1Clean(curRows,:);
    curM2 = m2Clean(curRows,:);
    
    W = learnOptimalOrthonormalTransformation(curM1, curM2);
    rowDelta = sqrt(sum((curM1 - curM2*W').^2,2));
    
    % Identify the keep rows for this iteration
    nKeepRows = nCleanRows - i;
    [~, sortOrder] = sort(rowDelta, 1, 'ascend');
    bestSortedRows = sortOrder(1:nKeepRows);
    curRows = sort(curRows(bestSortedRows));
    
end

alignRowsEnum = cleanRows(curRows)';
nAlignRows = length(alignRowsEnum); 

alignRows = false(size(m1,1), 1); 
alignRows(alignRowsEnum) = true; 

if nAlignRows < nStableRows
    warning('Too many small value rows: Unable to return the requested number of alignment rows.');
end
if nAlignRows < size(m1,2)
    warning('Number of alignment rows is less than the number of latent variables in the model.');

end