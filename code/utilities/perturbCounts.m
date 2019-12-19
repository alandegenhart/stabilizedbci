function trialData = perturbCounts(trialData, pertTrials, pertParams)
% Applies perturbations to raw counts recorded during an experiment to
% simulate instabilities. 
%
% Usage: trailData = perturbCounts(trialData, pertTrials, pertParams)
%
% Inputs:
%
%   trialData - a structure of trialData
%
%   pertTrials - indices of trials the perturbation should be applied to
%
%   pertParams - parameters specifying the perturbation that should be 
%   applied. See generatePertMatrices.m for more details. 
%
% Outputs: 
%
%   trialData - the trialData structure with the field 'pertBinCounts'
%   added.  If no perturbation has been applied to a trial, pertBinCounts
%   will just be equal to binCounts. However, if a pertubation has been
%   applied, then pertBinCounts will be the perturbed version of the 
%   original counts. 
%
% Author: William Bishop, bishopw@janelia.hhmi.org

% Generate matrices for applying the perturbation 
pertMats = generatePertMatrices(pertParams);
zeroAndPermute = pertMats.S*pertMats.P;
zeroAndOffset = pertMats.S*pertMats.O;

nTrials = length(trialData); 
for tI = 1:nTrials
    if any(tI == pertTrials)
        trialData(tI).pertBinCounts = bsxfun(@plus, zeroAndPermute*trialData(tI).binCounts', zeroAndOffset)';
    else
        trialData(tI).pertBinCounts = trialData(tI).binCounts; 
    end
end