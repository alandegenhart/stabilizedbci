function [C, psi, preComp] = fastFAMStep(blockEStr, blockXData, preComp)
% Computes the M-step of EM for fitting FA models. 
%
% Usage: fastFAMStep(blockEStr, blockXData, preComp)
%
% Inputs:
%
%   blockEStr: The output from the E-Step.  See fastFAEStep.m for more
%   details.
%
%   blockXData: The x data in block form.  Note that we assume x has been
%   mean centered. 
%
%   preComp: Saved values from an earlier call to this function. Saving
%   and reusing precomputed values can reduce computation.
%
% Outputs:
%
%   C, psi: The estimated values of C and psi.  See the function fitFA.m
%   for a description of FA model parameters.
%
%   preComp: pre-computed values which can be saved and used in later calls
%   to this function.
%
% Author: William Bishop, bishopw@janelia.hhmi.org

% Do precomputations, if needed 
if nargin < 3
    
    % We compute C, row-wise, in blocks; first we determine what our row-wise
    % blocks are
    xBlockInds = {blockXData.blockXInds};
    rowWiseBlockStr = findMinimalDisjointIntegerSets(xBlockInds, 'NO_ASSERT', true);

    nVarBlocks = length(rowWiseBlockStr); 
    for b = 1:nVarBlocks   
        blockCellInds = rowWiseBlockStr(b).cellInds;
        nCellInds = length(blockCellInds);
        
        rowWiseBlockStr(b).xBlockCols = cell(1, nCellInds); 
        for c = 1:nCellInds
            curCellInd = blockCellInds(c); 
            [~, rowWiseBlockStr(b).xBlockCols{c}] = ...
                intersect(blockXData(curCellInd).blockXInds, rowWiseBlockStr(b).vls);
        end
    end
    
    % Now we figure out the number of observed variables
    nVars = max([rowWiseBlockStr.vls]);
    
    % Now we figure out the number of latents
    nLatents = length(blockEStr(1).postCov);
    
    % Precompute square of observations
    obsSq =  nansum(blockStr2SmpMat(blockXData).^2)';
    
    % Precompute number of observations for each variable
    nObs = nansum(~isnan(blockStr2SmpMat(blockXData)))'; 
else
    rowWiseBlockStr = preComp.rowWiseBlockStr;
    nVars = preComp.nVars;
    nLatents = preComp.nLatents;
    obsSq = preComp.obsSq;
    nObs = preComp.nObs;
end

if nargout == 3
    preComp.rowWiseBlockStr = rowWiseBlockStr; 
    preComp.nVars = nVars; 
    preComp.nLatents = nLatents; 
    preComp.obsSq = obsSq; 
    preComp.nObs = nObs; 
end

C = nan(nVars, nLatents); 
psi = nan(nVars, 1);

% Now we compute C and and take care of some of the computations for psi in blocks, row-wise  
nVarBlocks = length(rowWiseBlockStr); 
for b = 1:nVarBlocks
    blockVarInds = rowWiseBlockStr(b).vls; 
    nBlockVars = length(blockVarInds); 
    
    xBlockInds = rowWiseBlockStr(b).cellInds; 
    nXBlocks = length(xBlockInds); 
    
    term1 = zeros(nBlockVars, nLatents); 
    term2 = zeros(nLatents, nLatents); 
    for xB = 1:nXBlocks
        
        curXBlockInd = xBlockInds(xB); 
        curXBlock = blockXData(curXBlockInd);
        curEStepBlock = blockEStr(curXBlockInd);
        
        curXBlockVarInds = rowWiseBlockStr(b).xBlockCols{xB};
        
        term1 = term1 + curXBlock.obsX(:, curXBlockVarInds)'*curEStepBlock.postMeans; 
        
        term2 = term2 + curXBlock.nSmps*curEStepBlock.postCov + ...
            curEStepBlock.postMeans'*curEStepBlock.postMeans;
    end
    
    C(blockVarInds,:) = term1/term2; 
    psi(blockVarInds) = obsSq(blockVarInds) - 1*sum(term1.*C(blockVarInds,:),2); 
end

psi = diag(psi./nObs);