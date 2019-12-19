function stabilizedFA = updateStabilizer(baseFA, y, varargin)
% This is a function to update a stabilizer using neural data.
%
% Usage: stabilizedFA = updateStabilizer(baseFA, y, varargin)
%
%   baseFA: The base stabilizer to update.  This should be a structure with
%           the fields d, C and psi as returned by fitBaseStabilizer.m
%
%   y: A cell of neural data to use in the update.  Each entry in the cell 
%   contains binned counts for one trial of shape n_neurons*n_bins.
%
% Optional Inputs: All optional inputs should be entered in string-value
% pair format.
%
%   N_FA_RESTARTS, MAX_N_ITS, LL_DIFF_THRESH, MIN_PRIV_VAR, C_INIT, 
%   PSI_INIT, VERBOSE: Optional inputs that are provided to fitFA.  
%   See fitFA.m for more details. Default values: 
%
%               N_FA_RESTARTS = 5;
%               MAX_N_ITS = 100000; 
%               LL_DIFF_THRESH = .00001; 
%               MIN_PRIV_VAR = .1;
%               C_INIT = [];
%               PSI_INIT = [];
%               VERBOSE = true;
%
%   ALIGN_N, ALIGN_TH: Optional inputs provided to alignLoadMatrices.  See
%   alignLoadingMatrices.m for more details.  Default values:
%
%               ALIGN_N = 0
%               ALIGN_TH = []
%
% Author: William Bishop, bishopw@janelia.hhmi.org

N_FA_RESTARTS = 5;
MAX_N_ITS = 100000; 
LL_DIFF_THRESH = .00001; 
MIN_PRIV_VAR = .1;
C_INIT = [];
PSI_INIT = [];
VERBOSE = true;

ALIGN_N = [];
ALIGN_TH = .01;

warnOpts(assignOpts(varargin)); 

if isempty(ALIGN_N)
    ALIGN_N = size(baseFA.C,1);
end

nLatents = size(baseFA.C,2); 
updateDataMat = [y{:}]';

faMdls = cell(1, N_FA_RESTARTS);
faLL = nan(1, N_FA_RESTARTS);
for mI = 1:N_FA_RESTARTS
            [fa.d, fa.C, fa.psi, ~, fa.diags] = ...
                fitFA(updateDataMat, nLatents,  'MAX_N_ITS', MAX_N_ITS, ...
                'LL_DIFF_THRESH', LL_DIFF_THRESH, ...
                'MIN_PRIV_VAR', MIN_PRIV_VAR, 'C_INIT', C_INIT, ...
                'PSI_INIT', PSI_INIT, 'VERBOSE', VERBOSE);
            faMdls{mI} = fa;
            faLL(mI) = fa.diags.ll(end);
end

% Pick the FA Model with the highest log-likelihood
% Pick the FA Model with the highest log-likelihood
[~, bestFAInd] = max(faLL);
stabilizedFA = faMdls{bestFAInd};

% Now we need to rotate the FA Model
[W, alignChs] = alignLoadingMatrices(baseFA.C, stabilizedFA.C, ALIGN_N, ALIGN_TH);
stabilizedFA.C = stabilizedFA.C*W';
stabilizedFA.alignChs = alignChs; 
