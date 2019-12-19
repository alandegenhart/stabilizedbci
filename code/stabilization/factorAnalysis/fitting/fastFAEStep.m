function blockEStr = fastFAEStep(blockStr, c, psi)
% Computes the E-step of EM for fitting FA models. 
%
% Usage: blockEStr = fastFAEStep(blockStr, c, psi)
%
% Inputs: 
%
%   blockXData: x data in block form.  Note that this function assumes
%   x data is already mean centered.
%
%   c, psi: loading matrix and matrix of private variances of the FA model,
%   respectively
%
% Outputs:
%
%   blockEStr: A structure with posterior means and posterior covariance
%   for each entry in blockStr.  Means will be contained in the field
%   'postMeans' and covariances in 'postCov'.
%
% Author: William Bishop, bishopw@janelia.hhmi.org
nLatents = size(c,2);
nBlocks = length(blockStr); 
for b = nBlocks:-1:1 % Work backwards to implicitly preallocate structure 
    
    curInds = blockStr(b).blockXInds;  
    
    blockC = c(curInds, :);
    blockPsi = psi(curInds, curInds); 
    
    coreComputation = blockC'/(blockC*blockC' + blockPsi);
    
    blockEStr(b).postMeans = blockStr(b).obsX*coreComputation';  
    blockEStr(b).postCov = eye(nLatents) - coreComputation*blockC;
end