function kf = kalmanFit(x, y, varargin)
% A function to learn the MLE parameters for a Kalman filter.  
%
% We consider a Kalman filter (KF) for the following dynamical system:
%
%   x_{t+1} = Ax_t + w_t
%       y_t = Cx_t + d + v_t ,
%
% where x_t is kinematics and y_t is neural data (which can be latent 
% state from a stabilizer), and we make the following distributional 
% assumptions:
%
%   x_0 ~ N(mu_1, V_1)
%   w_t ~ N(0, Q)
%   y_t ~ N(0, R)
%
% Given a set of training kinematics and neural data this function will
% form MLE estimates of the parameters of the dynamical system. 
%
% Usage: kf = kalmanFit(x, y)
%
% Inputs: 
%
%       x - a cell of length N.  Each entry contains a the latent state for a trial as a 
%       matrix of size xDim by T, where xDim is the dimensionality of x and T is the
%       number of steps in a particular trial
%
%       y - a cell of length N.  Each entry contains the observations for the corresponding 
%       trial in x as a matrix of size yDim by T, where yDim is the dimensionality of y.
%
% Optional Inputs: All optional inputs should be entered in string-value
% pair format. 
%   
%
%       FIT_STATE_NOISE - True if state noise covariance matrix, Q, should 
%       be fit; if false the user should specify STATE_NOISE_VAR, and the state 
%       noise covariance matrix, Q, will be assigned to be a diagonal matrix 
%       with the value STATE_NOISE_VAR along its diagonal.  Default: true
%
%       STATE_NOISE_VAR: The variance for each state variable if 
%       FIT_STATE_NOISE is false. Default: []
%
% Outputs: 
%
%   kf - a structure with the Kalman filter parameters.  It will contain
%   the fields d, A, C, mu_1, V_1, Q and R. 
% 
% Author: William Bishop, bishopw@janelia.hhmi.org
%

FIT_STATE_NOISE = true;
STATE_NOISE_VAR = [];
warnOpts(assignOpts(varargin)); 

% Remove empty trials
emptyTrials = cellfun('isempty', x); 
x(emptyTrials) = []; 
y(emptyTrials) = []; 

% Learn distribution of initial states
x0 = cellfun(@(xCell) xCell(:,1), x, 'UniformOutput', false); 
x0 = [x0{:}];
mu_1 = mean(x0,2); 
V_1 = cov(x0'); 

% Learn transition matrix
xWithTwoSteps = cellfun(@(xCell) size(xCell,2), x) > 1; 
xWithTwoStepsCell = x(xWithTwoSteps); 

xFirst = cellfun(@(xCell) xCell(:,1:end-1), xWithTwoStepsCell, 'UniformOutput', false); 
xSecond = cellfun(@(xCell) xCell(:,2:end), xWithTwoStepsCell, 'UniformOutput', false); 
xFirst = [xFirst{:}]; 
xSecond = [xSecond{:}]; 
A = xSecond/xFirst;

% Learn covariance matrix for state transition noise
if FIT_STATE_NOISE
    xSecondPred = A*xFirst; 
    xDelta = xSecond - xSecondPred; 
    Q = cov(xDelta');
else
    Q = diag(ones(size(xFirst,1),1)) * STATE_NOISE_VAR;
end

% Learning loading matrix & offset vector
allY = [y{:}]; 
allX = [x{:}]; 
allX =[allX; ones(1,size(allX,2))]; 
Cd = allY/allX; 
C = Cd(:,1:end-1); 
d = Cd(:,end);

% Learn the covariance matrix for the observation noise
yPred = Cd*allX; 
yDelta = allY - yPred; 
R = cov(yDelta'); 

% Return the KF structure
kf.d = d;
kf.A = A; 
kf.C = C;
kf.Q = Q; 
kf.R = R; 
kf.mu_1 = mu_1; 
kf.V_1 = V_1; 
