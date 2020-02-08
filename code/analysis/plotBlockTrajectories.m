% plotBlockTrajectories  Plot cursor trajectories for a single experimental
% block.
%
% Usage:
% plotBlockTrajectories(D)
%
% Inputs:
%   D    A single element of the structure array returned by
%        processTrialData
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

function plotBlockTrajectories(D)

% Get color info structure
C = defineColorInfo();

% Plot targets
r = D.TD(1).tgtRad;  % Target radius (same for all targets)
plotTargets(r, C)

% Iterate over all trials, plot each trajectory
for i = 1:D.nTrials
    % Get color
    [col, ~] = getColor(D.TD(i).tgtPos', C);
    
    % Get line style
    if logical(D.TD(i).success)
        lineStyle = '-';
    else
        lineStyle = '--';
    end
    
    % Plot trajectory
    plot(D.TD(i).refPos(:, 1), D.TD(i).refPos(:, 2), ...
        'color', col, ...
        'LineStyle', lineStyle) 
end

% Format plot
axLim = 150 * [-1, 1];
set(gca, ...
    'XLim', axLim, ...
    'YLim', axLim, ...
    'Box', 'on', ...
    'XTick', [], ...
    'YTick', [])

% Plot scale bar
scaleBar(gca, 50, 'mm', 'yScale', false, 'xScale', true, 'fontSize', 12)

end


function [colDark, colLight] = getColor(pos, C)
% Get color for associated target position.
%
% Usage:
% [colDark, colLight] = getColor(pos, C)
%

% Round target position to one decimal place
pos = round(pos, 1);

% Find matching indices for target position and get associated color
colIdx = ismember(C.targPos, pos, 'rows');
colDark = C.colDark(colIdx, :);
colLight = C.colLight(colIdx, :);

end


function plotTargets(r, C)
% Plot standard target configuration.
%
% Usage:
% plotTargets(r, C)
%

% Define target positions
targPos = getTargetPositions();

% Iterate over targets and plot
for i = 1:size(targPos, 1)
    [~, col] = getColor(targPos(i, :), C);
    rectangle( ...
        'Position', [targPos(i, 1) - r, targPos(i, 2) - r, r*2, r*2], ...
        'EdgeColor', 'k', ...
        'FaceColor', col, ...
        'Curvature', [1, 1])
end
end


function targPos = getTargetPositions()
% Get all possible locations of the targets. While normally this would be
% done by looking over the data directly, in this case we know exactly
% where the targets are so they can be defined explicitly.
%
% Usage:
% [targPos] = getTargetPositions()
%

r = 125;  % Target configuration radius (center-to-target distance)
ang = 0:45:315;
targPos = r * [cosd(ang)' sind(ang)'];
targPos = round(targPos, 1);  % Match precision of saved target positions

end


function C = defineColorInfo()
% Define mapping from target positions to color.
%
% Usage:
% [C] = defineColorInfo()
%

% Define color parameters
targPos = getTargetPositions();
nTarg = size(targPos, 1);
h = linspace(0, 1, nTarg+1)';
h = h(1:nTarg);
sDark = 1;
sLight = 0.4;
vLight = 1;
vDark = 0.75;

% Define HSV colors and convert to RGB
colHSVLight = [h, sLight * ones(nTarg, 1), vLight * ones(nTarg, 1)];
colHSVDark = [h, sDark * ones(nTarg, 1), vDark * ones(nTarg, 1)];
colLight = hsv2rgb(colHSVLight);
colDark = hsv2rgb(colHSVDark);

% Add info to output structure
C.targPos = targPos;
C.colLight = colLight;
C.colDark = colDark;

end