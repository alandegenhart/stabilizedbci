function [ci] = confidenceInterval(x,varargin)
% confidenceInterval        Calculate confidence interval for data
%
% [ci] = confidenceInterval(x)
%
% Plot confidence interval for the provided data.
%
%

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