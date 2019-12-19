function [ll, preComp] = faDataLogLikelihood(blockXData, c, psi, preComp)
% Calculates the likelihood of observed data under an FA model.
%
% Usage: [ll, preComp] = faDataLogLikelihood(blockXData, c, psi, preComp)
%
% Inputs: 
%
%   blockXData: x data in block form.  Note that this function assumes
%   x data is already mean centered.
%
%   c, psi: loading matrix and matrix of private variances of the FA model,
%   respectively
%
%   preComp: Optional inputs, pre-computations saved from an earlier call
%   to this function.
%
% Outputs:
%
%   ll: The log-likelihood
%
%   preComp: A structure of saved values which can be reused in future
%   calls to this function, potentially saving computational time
%
% Author: William Bishop, bishopw@janelia.hhmi.org
if nargin < 4
    preComp.nBlocks = length(blockXData); 
    for b = 1:preComp.nBlocks 
        nBlockVars = length(blockXData(b).blockXInds); 
        preComp.piConsts(b) = -.5*blockXData(b).nSmps*nBlockVars*log(2*pi);
        
        % For the way this function uses cov, we want to normalize by N and
        % not N - 1
        nSmps = blockXData(b).nSmps;
        if nSmps > 1
            preComp.covs{b} = ((nSmps-1)/nSmps)*cov(blockXData(b).obsX);
        else
            preComp.covs{b} = zeros(nBlockVars, nBlockVars); 
        end
        
    end
end

nBlocks = preComp.nBlocks;
piConsts = preComp.piConsts;
covs = preComp.covs;


ll = sum(piConsts); 

for b = 1:nBlocks
    
    nBlockSmps = blockXData(b).nSmps;
    blockVarInds = blockXData(b).blockXInds;
    curC = c(blockVarInds, :);
    curCov = curC*curC' + psi(blockVarInds, blockVarInds);

    ll = ll -.5*nBlockSmps*(log(det(curCov)) + trace(curCov\covs{b}));
end