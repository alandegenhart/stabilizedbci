% calculatePerformanceStatistics  Calculate evaluation block statistics.
%
% This function displays the statistics for the evaluation blocks of the
% provided session data.
%
% Usage:
%   calculatePerformanceStatistics(D)
%
% Inputs:
%   D   Block structure array (output from processTrialData)
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

function calculatePerformanceStatistics(D)

% Get evaluation blocks
D = D(logical([D.evaluationBlock]));

fprintf('\n---------------------------------------------------------------------\n')
fprintf('Block name \t\t SR (%%) \t AT (s) \t TAR (targ/s)\n')
fprintf('---------------------------------------------------------------------\n')
% Iterate over blocks
nBlocks = length(D);
for i = 1:nBlocks
    % Get success and acquisition time data
    success = logical(D(i).successCode);
    acquireTime = D(i).acquireTime;
    
    % Limit baseline washout block to 128 trials
    if strcmp(D(i).blockNotes, 'Baseline Washout')
        success = success(1:128);
        acquireTime = acquireTime(1:128);
    end
    
    % Calculate success rate, acquisition time, and target acquisition rate
    mean_sr = mean(success) * 100;  % Convert to percentage
    mean_at = mean(acquireTime(success))/1000;  % Convert to seconds
    tar = sum(success)/(sum(acquireTime)/1000);  % targets/second
    
    % If success rate is less than 50%, set acquisition time to be NaN
    if mean_sr < 50
        mean_at = nan;
    end
    
    nameStr = getShortBlockName(D(i).blockNotes);
    fprintf('%s \t %6.2f \t %0.2f \t\t %0.2f\n', ...
        nameStr, mean_sr, mean_at, tar)
end

fprintf('---------------------------------------------------------------------\n\n')
end


function blockName = getShortBlockName(blockName)
% getShortBlockName  Get shortened block name for displaying results.
%
% Usage:
% [blockName] = getShortBlockName(blockName)
%

% Note: don't need to change baseline washout
switch blockName
    case 'Baseline Evaluation'
        blockName = 'Baseline Eval.   ';
    case 'Stitching Evaluation'
        blockName = 'Stabilizer Eval. ';
    case 'Perturbation Evaluation'
        blockName = 'Instability Eval.';
    case 'Baseline Washout'
        blockName = 'Baseline Washout ';
end
end