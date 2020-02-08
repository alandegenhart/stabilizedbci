% This script replicates the main analysis used to create Figure 4 of the
% manuscript.
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

%% ========================================================================
% Add the required paths to use the stabilization code
startDir = pwd;
cd('code'); 
addStabilizationProjectPaths();
cd(startDir); 

%%  =======================================================================
% Load neural data for the example session provided.

dataFile = fullfile('data', '20160325');
data = load(dataFile); 

%% ========================================================================
% Create plots summarizing performance

% Process trial data and print session summary statistics
D = processTrialData(data.trialData);
calculatePerformanceStatistics(D)

% Plot session performance (Fig. 4b)
[F1] = plotSessionPerformance(D, 'Example session - L20160325');
% Plot cursor trajectories (Fig. 4c)
[F2] = plotCursorTrajectories(D, 'Example session trajectories - L20160325');
