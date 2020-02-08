% plotCursorTrajectories  Plot cursor trajectories for stabilization
% session.
%
% Usage:
% [F] = plotCursorTrajectories(D)
%
% Inputs:
%   D         Structure array
%   titleStr  Title of figure
%
% Outputs:
%   F         Figure handle
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

function [F] = plotCursorTrajectories(D, titleStr)

% Get blocks to plot.  Find the indices for the evaluation blocks
% (baseline, instability, and stabilization) and get the first 3 stabilizer
% blocks.
blockNames = {D.blockNotes};
baselineEvalInd = find(strcmp(blockNames, 'Baseline Evaluation'));
stabilizerInd = find(strcmp(blockNames, 'Stitching'));
stabilizerEvalInd = find(strcmp(blockNames, 'Stitching Evaluation'));
instabilityEvalInd = find(strcmp(blockNames, 'Perturbation Evaluation'));

% Update block names. These are used as the title for each panel.
D(stabilizerInd(1)).blockNotes = 'Instability';
for i = 2:length(stabilizerInd)
    D(stabilizerInd(i)).blockNotes = sprintf('Stabilizer Update %d', i-1);
end
D(stabilizerEvalInd).blockNotes = 'Stabilizer Evaluation';
D(instabilityEvalInd).blockNotes = 'Instability Evaluation';
D = D([ ...
    baselineEvalInd, ...
    stabilizerInd(1:3), ...
    stabilizerEvalInd, ...
    instabilityEvalInd]);
nBlocks = length(D);

% =========================================================================
% Setup figure
nCol = nBlocks;
nRow = 1;
axW = 150;
axH = 150;
axSp = 20;
xMargin = [75, 75];
[fW,fH,Ax] = calcFigureSize(nRow, nCol, axW, axH, axSp, ...
    'xMargin', xMargin);
F = figure('Position',[100 100 fW fH]);
set(F, 'color', 'w')

% Plot figure title
plotTitle(titleStr)

% =========================================================================
% Plot trajectories

for i = 1:nBlocks
    subplotSimple(nRow, nCol, i, 'Ax', Ax); hold on;
    plotBlockTrajectories(D(i))
    title(D(i).blockNotes, 'FontSize', 14)
end