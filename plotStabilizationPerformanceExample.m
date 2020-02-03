% This script replicates the main analysis used to create Figure 4 of the
% manuscript.

%%  =======================================================================
% Add the required paths to use the stabilization code
startDir = pwd;
cd('code'); 
addStabilizationProjectPaths();
cd(startDir); 

%%  =======================================================================
% Load neural data for the example session provided.

dataFile = fullfile('data', '20160325');
data = load(dataFile); 

%% =======================================================================
% Create plots summarizing performance

% Process trial data
processTrialData(data.trialData)

% Analysis:
% - Main function to calculate success rate/target acquisition rate as a
% function of target block
% - Summarize overal stabilization performance for the session

% Plots to generate:
% - Performance vs block (Fig. 4b)
% - Trajectories for individual blocks (Fig. 4c)
