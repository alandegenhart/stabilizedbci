function setStr = findMinimalDisjointIntegerSets(inputCell, varargin)
% Given sets of possibly overlapping integers, this function will return another
% collection of sets so that intersection of any two of the returned sets
% is the null set and the union of the returned sets equals the union of
% the original sets. This collection will be the smallest collection 
% possible with this property. 
%
% Note: This function is optimized for speed and only operates on integers.
%
% Usage: setStr = findMinimalDisjointIntegerSets(inputCell, varargin)
%
% Inputs: 
%
%   inputCell - A cell array, where each entry in the array contains a
%   vector of integers.  
%
% Optional Inputs: All optional inputs should be entered in string value
% pair format. 
%
%   NO_ASSERT - If true, assertions will be skipped.  Good for speed but
%   may not catch unexpected input.  Default: false
%
% Outputs: 
%
%   setStr - a structure of length S, where S is the numbered of returned
%   sets.  It will have the following fields:
%
%       vls - the values in the set
%
%       cellInds - the indices into inputCell of the original sets that
%       contained these values. 
%
% Author: William Bishop, bishopw@janelia.hhmi.org

NO_ASSERT = false; 
warnOpts(assignOpts(varargin)); 

if ~NO_ASSERT
    assert(all(cellfun(@(x) all(isRoundNum(x)), inputCell)), 'Entries of inputCell must contain integers.'); 
    assert(all(cellfun(@(x) isvector(x), inputCell)), 'Entries of inputCell must be vectors.'); 
end

cellLength = length(inputCell); 

minVl = min(cellfun(@(x) min(x), inputCell));
maxVl = max(cellfun(@(x) max(x), inputCell));

uniqueVls = minVl:1:maxVl;
nUniqueVls = maxVl - minVl+1;

indMatrix = false(nUniqueVls, cellLength); 
for i = 1:cellLength
    indMatrix(inputCell{i}-minVl+1,i) = true; 
end

badRows = sum(indMatrix,2) == 0;
uniqueVls(badRows) = [];
indMatrix(badRows,:) = []; 

[uniqueRows, ~, rowMap] = unique(indMatrix, 'rows'); 
nUniqueRows = size(uniqueRows,1); 

for r = nUniqueRows:-1:1
    setStr(r).vls = uniqueVls(rowMap == r);
    setStr(r).cellInds = find(uniqueRows(r,:)); 
end

% Handle the case of empty input
if nUniqueRows == 0
    setStr.vls = []; 
    setStr.cellInds = []; 
end

