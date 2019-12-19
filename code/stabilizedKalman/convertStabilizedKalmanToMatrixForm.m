function mMatrices = convertStabilizedKalmanToMatrixForm(dc, varargin)
% This is a function to take a stabilized Kalman decoder in standard form
% (that is in the form returned by stabilizedKalmanFit) and put it into 
% "M-matrix" form so that we can decode with:
%
%     x_t = M_2y_t + M_1x_{t-1} + M_0, 
%
% where y_t is observed data and x_t is decoded state.  For the first step,
% x_{t-1} will be set to mu_1.  We are able to put the stabilizer and
% Kalman filter into this matrix form by assuming that the Kalman gain
% has converged. As such, this functin will compute the converged Kalman
% gain. 
%
% Note: 
%
% Usage: mMatrices = convertStabilizedKalmanToMatrixForm(dc)
%
% Inputs: 
%
%   dc - a structure with the model parameters for the KalmanFA decoder.
%   See stabilizedKalmanFit.m for more information. 
%
% Optional Inputs:
%
%   COMPUTE_METHOD: As noted in the publication, this function originally 
%   omitted  multiplying by the state transition matrix when calculating M_1.
%   Additionally, the iterations for calcuating the converged Kalman gain 
%   were initialized incorrectly in the original training code.  Both of
%   these have been corrected in the default version of the code. However,
%   to obtain the M-matrices computed in the closed-loop experiments, set
%   'COMPUTE_METHOD' to 'OLD'. Default: 'CORRECT'
%
%   MAX_GAIN_STEPS: The max number of iterations to allow for calculating
%   the converged Kalman gain. Default: 100000
%
%   STOP_TOL: The stopping tolerance when computing the converged Kalman
%   gain. Default 1E-8
%
% Outputs: 
%
%   mMatrices - a structure with the mu_1, M0, M1 and M2 matrices as fields
%
% Author: William Bishop, bishopw@janelia.hhmi.org

COMPUTE_METHOD = 'CORRECT'; % 'CORRECT' or 'OLD'
MAX_GAIN_STEPS = 100000;
STOP_TOL = 1E-8;
    
warnOpts(assignOpts(varargin));

% Calculate beta
xDim = size(dc.kf.A,1);
beta = dc.fa.C'/(dc.fa.C*dc.fa.C' + dc.fa.psi);

if strcmp(COMPUTE_METHOD, 'CORRECT')
    % Calculate the converged Kalman gain, correcting for how we initialize
    C = dc.kf.C;
    R = dc.kf.R;
    A = dc.kf.A;
    Q = dc.kf.Q;
    V1 = dc.kf.V_1;
    
    curStep = 1;
    prevGain = inf(size(C'));
    delta = inf;
    oneStepV = V1;
    while curStep <= MAX_GAIN_STEPS && delta > STOP_TOL
        curGain = oneStepV*C'/(R+C*oneStepV*C');
        
        curV = oneStepV - curGain*C*oneStepV;
        oneStepV = A*curV*A' + Q;
        
        delta = sum(sum((curGain - prevGain).^2));
        prevGain = curGain;
    end
    if delta >= STOP_TOL
        warning('Unable to compute converged Kalman gain.');
    end
    dc.kf.cGain = curGain;
    
    % Calculate correct m-matrices, correcting for how we initialize
    mMatrices.M2 = dc.kf.cGain*beta;
    mMatrices.M1 = (eye(xDim) - dc.kf.cGain*dc.kf.C)*dc.kf.A;
    mMatrices.M0 = -dc.kf.cGain*(dc.kf.d + beta*dc.fa.d);
    mMatrices.mu_1 = dc.kf.mu_1;
    
elseif strcmp(COMPUTE_METHOD, 'OLD')
    C = dc.kf.C;
    R = dc.kf.R;
    A = dc.kf.A;
    Q = dc.kf.Q;
    V_1 = dc.kf.V_1;

    curStep = 1;
    prevV = V_1;
    prevGain = inf(size(C')); 
    delta = inf; 
    while curStep <= MAX_GAIN_STEPS && delta > STOP_TOL
        curStep = curStep + 1; 
        oneStepV = A*prevV*A' + Q;
        curGain = oneStepV*C'/(R+C*oneStepV*C');
        prevV = oneStepV - curGain*C*oneStepV;
        delta = sum(sum((curGain - prevGain).^2));
        prevGain = curGain;
    end
    if delta >= STOP_TOL
        warning('Unable to compute converged Kalman gain.'); 
    end
    dc.kf.cGain = curGain;
    
    mMatrices.M2 = dc.kf.cGain*beta;
    mMatrices.M1 = (eye(xDim) - dc.kf.cGain*dc.kf.C);
    mMatrices.M0 = -dc.kf.cGain*(dc.kf.d + beta*dc.fa.d);
    mMatrices.mu_1 = dc.kf.mu_1;
else
    error(['The COMPUTE_METHOD ', COMPUTE_METHOD, ' it not recognized.']);
end






