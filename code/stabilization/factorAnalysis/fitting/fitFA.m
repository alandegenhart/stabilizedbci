function [d, c, psi, conv, diags] = fitFA(x, nLatents, varargin)
% This is a function to fit an FA model. 
%
% In particular, we model:
%
%   x_t = C * l_t + d + ep_t, 
%
% where we assume:
%
%   l_t ~ N(0, I)
%
%   ep_t ~ N(O, psi), where psi is a diagonal matrix.  Given training
%   data, x, in a matrix this function will fit C, d and psi using
%   expectation maximization.  This function can also handle missing
%   values in x, but this functionality is not used for stabilization. 
%
% Usage: [d, c, psi, conv, diags] = fitFA(x, nLatents, varargin)
%
% Inputs: 
%
%   x - an n by p data matrix.  Each row is a sample. Missing values are
%   indicated by nan. 
%
%   nLatents - the number of latents to fit the model for. 
%
% Optional Inputs: All optional inputs should be entered in string-value
% pair format.
%
%   MAX_N_ITS - The max number of EM iterations that the method should use
%   when fitting the model via EM. Default: 100000. 
%
%   LL_DIFF_THRESH - The stopping criterion for fitting the model via EM.
%   Default: 1E-8.
%
%   MIN_PRIV_VAR - After the M-step of each fitting iteration, we can
%   enfore that the private variance values (the diagonal values along psi)
%   are all above a certain threshold.  MIN_PRIV_VAR gives this threshold. 
%   Default: .01
%
%   C_INIT - If provided, the initial value of C to start EM with.  If
%   empty, a random C will be used.  Default: []
%
%   PSI_INIT - If provided, the initial value of psi to start EM with.  If
%   empty, a random psi will be used.  Default: []
%
%   VERBOSE - True if this function should print periodic status updates to
%   screen. Default: true
%
% Outputs: 
%
%   d - A P by 1 vector of the estimated model mean. 
%
%   c - A P by L matrix giving the estimated loading matrix
%
%   psi - A P by P diagonal matrix of the estimated private noise
%   variances. 
%
%   conv - True if the final EM algorithm converged
%
%   diags - A structure with diagnostic information.  It will contain the
%   following: 
%
%       cInit, psiInit - the initial values of the c and psi parameters
%
%       ll - The log-likelihhod after each EM iteration. ll(i+1) gives the
%       log-likelihood after the i^th M step.  ll(1) gives the starting 
%       log-likelihood with the initial parameters. 
%
%       fitFcn - The name of the fitting function. In this case,
%       'fitFAWithMissingVls'
%
% Author: William Bishop, bishopw@janelia.hhmi.org

MAX_N_ITS = 100000;
LL_DIFF_THRESH = 1E-8;
VERBOSE = true;
MIN_PRIV_VAR = .01; % Used to combat Heywood cases
C_INIT = [];
PSI_INIT = [];

warnOpts(assignOpts(varargin));

% Setup the function we use for our stopping criteria
stopFcn = @(cL, pL, iL)  (cL - pL);

% See if diags were requested
diagsReq = nargout == 5;

% See how many observed variables we are dealing with
nObsVars = size(x,2);

% Calculate our mean
d = nanmean(x)';

if nLatents == 0 % Handle easy case
    c = zeros(nObsVars,1);
    
    psi = nanvar(x);
    psi(psi < MIN_PRIV_VAR) = MIN_PRIV_VAR;
    psi = diag(psi);
    
    diags = struct([]); 
    
    conv = true;
else % Run EM if there are any latents
    
    % Subtract the mean from x
    xCtr = bsxfun(@minus, x, d');
    
    % Arrange centered data into blocks according to what variables are
    % observed. This is for computational efficiency. 
    xCtrBlocked = smpMat2BlockStr(xCtr);

    % Initialize estimates for c and psi
    if isempty(C_INIT) || isempty(PSI_INIT)
            if VERBOSE
                disp('Initializing randomly.'); 
            end
            % Initialize c matrix 
            cInit = randn(nObsVars, nLatents);
            % Initialize Psi 
            psiInit = diag(nanvar(x)) - diag(diag(cInit*cInit'));
    end
    
    % Override initial parameters if the user has specified these
    if ~isempty(C_INIT)
        cInit = C_INIT;
    end
    if ~isempty(PSI_INIT)
        psiInit = PSI_INIT;
    end
    
    if diagsReq
        diags.cInit = cInit; 
        diags.psiInit = psiInit; 
    end
     
    % Give user some more feedback
    if VERBOSE
        disp(['Done with initialization.  Fitting with EM.']); 
    end
    
    % Allocate space for diagnostics if requested
    if diagsReq
        diags.ll = nan(1, MAX_N_ITS+1);
    end
    
    c = cInit;
    psi = psiInit;
    
    % Enforce private noise floor on psi
    psiDiag = diag(psi);
    psiDiag(psiDiag < MIN_PRIV_VAR) = MIN_PRIV_VAR;
    psi = diag(psiDiag);
    
    % Calculate initial log-likelihood
    [initLL, llPreComp] = faDataLogLikelihood(xCtrBlocked, c, psi);

    if diagsReq
        diags.ll(1) = initLL;
    end
    prevLL = initLL;
    
    curIt = 1;
    stopCrit = inf;
    while curIt <= MAX_N_ITS && stopCrit > LL_DIFF_THRESH
        
        % Run the E-Step
        blockEStr = fastFAEStep(xCtrBlocked, c, psi);
        
        % Run the M-Step
        if curIt == 1
            [c, psi, preComp] = fastFAMStep(blockEStr, xCtrBlocked);
        else
            [c, psi] = fastFAMStep(blockEStr, xCtrBlocked, preComp);
        end
        
        % Enforce private noise floor 
        psiDiag = diag(psi);
        psiDiag(psiDiag < MIN_PRIV_VAR) = MIN_PRIV_VAR;
        psi = diag(psiDiag);
        
        curLL = faDataLogLikelihood(xCtrBlocked, c, psi, llPreComp);
        diags.ll(curIt + 1) = curLL;
    
        if curLL == -inf
            stopCrit = inf;
        else
            stopCrit = feval(stopFcn, curLL, prevLL, initLL); 
        end 
        prevLL = curLL;
        
        % Display Results
        if VERBOSE && mod(curIt, 100) == 0
            fprintf('EM Iteration: %5i LL: %0.8g Stop Crit: %0.12g \r', curIt, curLL, stopCrit);
        end
        
        % Iterate counter
        curIt = curIt + 1;
    end
    
    % Report if we converged or not
    conv = stopCrit <= LL_DIFF_THRESH;
    diags.ll(curIt+1:end) = [];
end

if diagsReq
    diags.fitFcn = 'fitFA'; 
end