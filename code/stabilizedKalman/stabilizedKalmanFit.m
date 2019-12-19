function mP = stabilizedKalmanFit(x, y, nLatents, varargin)
% A function to fit an initial stabilizer and Kalman filter to 
% calibration data.  
%
% Usage: mP = stabilizedKalmanFit(x, y, nLatents, varargin)
%
% Inputs:
%
%       x - trainining kinematic data.  This should be a cell of length N
%       where N is the number of training trials. Each entry contains
%       kinematic data of size xDim by T, where xDim is the number of 
%       kinematic variables and T is the number of steps in a particular 
%       trial.
%
%       y - training neural data.  This should be a cell of length N.  Each 
%       entry contains the observations for the corresponding trial in x as
%       a  matrix of size yDim by T, where yDim is the dimensionality of y.
%
% Optional Inputs: All optional inputs should be entered in string-value
% pair format.
%
%       N_FA_RESTARTS, MAX_N_ITS, LL_DIFF_THRESH, MIN_PRIV_VAR, C_INIT,
%       PSI_INIT, VERBOSE: These are optional inputs when fitting the base 
%       stabilizer.  See the function fitBaseStabilizer.m for more information. 
%       Default values:
%
%           N_FA_RESTARTS = 5;
%           MAX_N_ITS = 100000;
%           LL_DIFF_THRESH = .00001; 
%           MIN_PRIV_VAR = .01; 
%           C_INIT = [];
%           PSI_INIT = [];
%           VERBOSE = true;
%
%
% Outputs:
%
%   mP - a structure with the model parameters for the KalmanFA decoder. It
%   will have the fields kf (containing the parameters for the kalman
%   filter, see the function kalmanFit for a description of the fields in
%   kf) and fa (containing parameters of the baseline stabilizer, see
%   fitFAWithMissingVls.m for more information). 
%
% Author: William Bishop, bishopw@janelia.hhmi.org

N_FA_RESTARTS = 5;
MAX_N_ITS = 100000;
LL_DIFF_THRESH = .00001;
MIN_PRIV_VAR = .01; 
C_INIT = [];
PSI_INIT = [];
VERBOSE = true;

FIT_STATE_NOISE = true;
STATE_NOISE_VAR = [];

warnOpts(assignOpts(varargin)); 

% Fit the baseline stabilizer
fa = fitBaseStabilizer(y, nLatents, 'N_FA_RESTARTS', N_FA_RESTARTS, ...
    'MAX_N_ITS', MAX_N_ITS, 'LL_DIFF_THRESH', LL_DIFF_THRESH, ...
    'MIN_PRIV_VAR', MIN_PRIV_VAR, 'C_INIT', C_INIT, 'PSI_INIT', PSI_INIT, ...
    'VERBOSE', VERBOSE); 

% Infer latents using the FA model
nTrials = length(y);
l = cell(1, nTrials);
for tI = 1:nTrials
    if ~isempty(y{tI})
        curL= inferFALatents(fa, y{tI}');
        l{tI} = curL';
    end
end

% Train the Kalman filter portion of the model
kf = kalmanFit(x, l, 'FIT_STATE_NOISE', FIT_STATE_NOISE, ...
    'STATE_NOISE_VAR', STATE_NOISE_VAR);

mP.fa = fa;
mP.kf = kf;
