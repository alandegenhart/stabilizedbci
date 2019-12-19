function xHat = stabilizedKalmanDecodeMatrixForm(mMatrices, data, varargin)
% A function to decode data with a stabilized Kalman decoder in matrix form.
%
% Usage: xHat = stabilizedKalmanDecodeMatrixForm(mMatrices, data, varargin)
%
% Inputs:
%
%   mMatrices - a structure with the M-matrices for the decoder.  See the
%   fuction convertStabilizedKalmanToMatriForm.m for more information. 
%
%   data - a cell of length N, where N is the number of trials.  data{n}
%   contains neural data, Y, for the n^th trial.  Y should be size yDim by T,
%   where yDim is the number of neurons and T is the number of steps
%   in the trial.
%
% Optional Inputs: All optional inputs should be given in string-value pair
% format.
%
%   VERBOSE - True if progress updates should be output to screen.
%             Default: true 
%
% Outputs:
%
%   xHat - a cell of length N containing the decoded state for each trial.
%   Each entry will contain a matrix of size x_dim by T, were xDim is the
%   dimensionality of the state (e.g., kinematics). 
%
% Author: William Bishop, bishopw@janelia.hhmi.org
%
VERBOSE = true;
warnOpts(assignOpts(varargin));

xDim = length(mMatrices.mu_1);
nTrials = length(data);

xHat = cell(1, nTrials);
for tI = 1:nTrials
    curData = data{tI};
    nSteps = size(curData,2);
    
    curXHat = nan(xDim,nSteps);
    prevX = mMatrices.mu_1;
    for sI = 1:nSteps
        stepX = mMatrices.M2*curData(:,sI) + mMatrices.M1*prevX + mMatrices.M0;
        curXHat(:, sI) = stepX;
        prevX = stepX;
    end
    xHat{tI} = curXHat;
    
    if VERBOSE && mod(tI,100) == 0
        disp(['Done decoding trial ', num2str(tI), ' of ', num2str(nTrials)]);
    end
end
