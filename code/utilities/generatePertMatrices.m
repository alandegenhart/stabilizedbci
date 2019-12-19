function pM = generatePertMatrices(pertParams)
% Given the values in a structure of perturbation parameters, this function
% will generate matrices S, P and O, so that if
%
%   y' = S*P*y + S*O, 
%
% then y' will be perturbed neural data generated from y according the
% values given in structure of perturbation parametres. In particular, the
% matrix P will switch (permute) channels, S silence channels and O will
% offset channels.
%
% Usage: pM = generatePertMatrices(pertParams)
%
% Inputs: 
%
%   pertParams - a structure of perturbation parameters.  It must have the
%   following fields: 
%
%       nArrayChs - the number of channels in y (e.g., the number of channels 
%       on the array) 
%
%       baseElectrodes - an array stating which electrodes are used by the
%       decoder before we simulate applying any perturbations. 
%
%       swapElectrodesMap - an array of size B by 2, where B is the number
%       of base electrodes.  The first column should be a copy of
%       baseElectrodes. Each row contains an electrode "swap" in the
%       convention that counts for electrodes in the second column will
%       be substituted for counts of the electrodes in the first column (so
%       to indicate an electrode should not be swapped, the entries in the
%       first and second column for that electrode's row should be the
%       same). We can swap in a heldout unit by using an electrode in the
%       second column that is not among the base electrodes. 
%       Note: It is not a true "swap" in the sense that counts for
%       electrodes in the first column *ARE NOT* in turn substituted for
%       counts of electrodes in the second. 
%
%       silentChs - a logical array of length B.  If silentChs(b) is true,
%       it indicates that the b^th base electrode (electrode indicated by the
%       b^th entry of baseElectrodes should be silenced).  Note (see
%       formula above) that silencing is applied after swapping and
%       offsets. 
%
%       offsets - an array of length B.  offsets(b) contains the offsets to
%       be applied to the counts for baseUnits(b).  Note that offsets are
%       applied to swapped counts and before silencing (see formula above).
%
% Outputs: 
%
%   pM - a structure with the fields S, P and O holding the above matrices.
%   
%
% Author: William Bishop, bishopw@janelia.hhmi.org

P = eye(pertParams.nArrayChs); 
for cI = 1:size(pertParams.swapElectrodesMap,1)
    P(pertParams.swapElectrodesMap(cI,1), pertParams.swapElectrodesMap(cI,1)) = 0; 
    P(pertParams.swapElectrodesMap(cI,1), pertParams.swapElectrodesMap(cI,2)) = 1; 
end

sDiag = ones(pertParams.nArrayChs,1); 
sDiag(pertParams.baseElectrodes) = double(~pertParams.silentChs); 
S = diag(sDiag); 

O = zeros(pertParams.nArrayChs,1); 
O(pertParams.baseElectrodes) = pertParams.offsets;

pM.P = P;
pM.S = S;
pM.O = O; 


