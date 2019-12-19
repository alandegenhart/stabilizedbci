function fa = fitBaseStabilizer(y, nLatents, varargin)
% Fits a base stabilizer to neural data.
%
% A base stabilizer is an FA model fit to data, so this function is a 
% wrapper for fitting FA models, which performs multiple random restarts
% and picks the FA model with the highest likelihood on the provided
% neural data. 
%
% Inputs:
%
%   y - a cell of neural data.  Each entry in the cell contains binned
%   counts for one trial of shape n_neurons*n_bins
%
%   nLatents - The number of latents in the fit FA models.
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
% Outputs:
%
%   fa - the fit FA model.  With the parameters d (mean vector), 
%   C (loading matrix) and psi (private variances).
%
% Author: wbishop@janelia.hhmi.org
%

N_FA_RESTARTS = 5;
MAX_N_ITS = 100000;
LL_DIFF_THRESH = .00001;
MIN_PRIV_VAR = .1; 
C_INIT = [];
PSI_INIT = [];
VERBOSE = true;

warnOpts(assignOpts(varargin));

yMat = [y{:}];

% Fit the FA Model, using random restarts
faMdls = cell(1, N_FA_RESTARTS);
faLL = nan(1, N_FA_RESTARTS);
for mI = 1:N_FA_RESTARTS
    [fa.d, fa.C, fa.psi, ~, fa.diags] = ...
        fitFA(yMat', nLatents, 'MAX_N_ITS', MAX_N_ITS, ...
        'LL_DIFF_THRESH', LL_DIFF_THRESH, 'MIN_PRIV_VAR', MIN_PRIV_VAR, ...
        'C_INIT', C_INIT, 'PSI_INIT', PSI_INIT, 'VERBOSE', VERBOSE);
    faMdls{mI} = fa;
    faLL(mI) = fa.diags.ll(end);
end

% Pick the FA Model with the highest log-likelihood
[~, bestFAInd] = max(faLL);
fa = faMdls{bestFAInd};