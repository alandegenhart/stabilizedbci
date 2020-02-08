% confidenceInterval  Calculate confidence interval for data.
%
% This function calculates the confidence interval for the provided data.
% Confidence intervals can be estimated using a t-distribution or by
% resampling. When resampling is specified, random samples are drawn with
% replacement from the input data.
%
% Usage:
%   [ci] = confidenceInterval(x)
%
% Inputs:
%   x       Input data array
%
% Optional inputs:
%   alpha       Alpha value used to calculate the confidence interval
%   method      Method used to calculate the confidence interval
%               ('bootstrap' or 'tDist')
%   statistic   String specifying the function name of the statistic to
%               calculate the confidence for
%   nRep        Number of repetitions to use (bootstrap method only)
%
% Outputs:
%   ci          Confidence interval
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

function [ci] = confidenceInterval(x,varargin)

% Parse optional arguments
alpha = 0.05;
method = 'bootstrap';
statistic = 'mean';
nRep = 10000;

assignOpts(varargin);

% Reshape to 1 x samples and remove NaNs
x = reshape(x,1,[]);
x = x(~isnan(x));
nSamp = length(x);

switch method
    case 'tDist'
        nDoF = nSamp - 1;

        % Calculate standard deviation and standard error
        stDev = std(x);
        stErr = stDev/sqrt(nSamp);

        % Calculate t-statistic
        int = [0 1] + ([1 -1] * (alpha/2));
        tStat = tinv(int,nDoF);

        % Calculate confidence interval
        ci = tStat * stErr;
        
    case 'bootstrap'
        fn = str2func(statistic);
        % Loop over iterations
        meanBS = nan(nRep,1);
        for i = 1:nRep
            % Draw nSamp samples with replacement and calculate mean
            ind = randi(nSamp,1,nSamp);
            meanBS(i) = fn(x(ind));
        end
        
        % Sort by asending statistic and calculate upper and lower bound
        meanBS = sort(meanBS);
        nPtsTail = round(nRep*alpha/2);
        lCI = meanBS(nPtsTail);
        uCI = meanBS(nRep - nPtsTail + 1);
        ci = [lCI uCI];
end