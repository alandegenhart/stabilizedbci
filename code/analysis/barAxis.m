% barAxis  Generate bracket-style axis bars.
%
% This function adds 'bracket-style' axes to the current axis. Multiple
% sets of brackets can be added to the same axis, effectively allowing the
% axis to be "broken."
%
% Usage:
%   barAxis(axType, axLim, axTick)
%
% Inputs:
%   axType      Axis to add bracket axis to ('x' or 'y')
%   axLim       Axis limits for each set of axis brackets to plot. Each set
%               set of axis brackets should be a separate row.
%   axTick      Tick marks for each axis
%
% Optional inputs:
%   axisTickLabels  Cell array of axis tick labels to use
%   plotTickLabels  Specifies whether or not to plot tick labels
%   tickSize        Length of tick marks
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

function barAxis(axType, axLim, axTick, varargin)

% Optional arguments
axisTickLabels = [];
plotTickLabels = true;
tickSize = 5;

assignOpts(varargin, false);

% If axis tick label strings are not specified, convert axis tick numbers
% to axis label strings.
nTick = size(axTick,2);
nAx = size(axLim,1);
if isempty(axisTickLabels)
    axisTickLabels = cell(nAx,nTick);
    for i = 1:nTick
        for j = 1:nAx
            axisTickLabels{j,i} = num2str(axTick(j,i));
        end
    end
end

% If axis tick labels is not a cell array of strings, convert appropriately
if isnumeric(axisTickLabels)
    tempAxisTickLabels = axisTickLabels;
    axisTickLabels = cell(nAx,nTick);
    for i = 1:nTick
        for j = 1:nAx
            axisTickLabels{j,i} = num2str(tempAxisTickLabels(j,i));
        end
    end
end

% Check to ensure specified tick labels are the same size as the ticks
if size(axTick,2) ~= size(axisTickLabels,2)
    error('Number of axis labels does not match the number of specified tick marks.')
end

% Get axis limits.  These will be used to determine the size of the tick
% marks
xLim = get(gca,'XLim');
yLim = get(gca,'YLim');

% Determine tick sizes.  These are be independent of the axis scale, and
% are a fixed percentage of the figure size.
rngW = xLim(2) - xLim(1);
rngH = yLim(2) - yLim(1);
axPos = get(gca,'Position');
axW = axPos(3);
axH = axPos(4);
figPos = get(gcf,'Position');
figW = figPos(3);
figH = figPos(4);
yTickSz = (tickSize * rngW) / (figW*axW - tickSize);
xTickSz = (tickSize * rngH) / (figH*axH - tickSize);

% Loop over axis label sets.  This allows multiple bar plots to be
% generated.
for axNo = 1:nAx
    % Determine X and Y axis data to plot depending on 'axType'
    tempAxLim = axLim(axNo,:);
    tempAxTick = axTick(axNo,:);
    switch lower(axType)
        case 'x'
            y = yLim(1) * ones(1,2);
            x = tempAxLim;
            xTick = [tempAxTick;tempAxTick];
            yTick = [y(1) y(1) - xTickSz]';
            
            % Update axis limits
            xLimNew = xLim;
            yLimNew = [(yLim(1) - xTickSz) yLim(2)];

            % Specify location of labels
            xTextLoc = tempAxTick;
            yTextLoc = (yLimNew(1) - xTickSz/2) * ones(1,nTick);

            % Set axis name to turn off and text alignment string
            axName = 'XAxis';
            horizAlignStr = 'center';
            vertAlignStr = 'top';
        case 'y'
            y = tempAxLim;
            x = xLim(1) * ones(1,2);
            xTick = [x(1) x(1) - yTickSz]';
            yTick = [tempAxTick;tempAxTick];
            
            % Update axis limits
            xLimNew = [(xLim(1) - yTickSz) xLim(2)];
            yLimNew = yLim;

            % Specify location of tick labels
            xTextLoc = (xLimNew(1) - yTickSz/2) * ones(1,nTick);
            yTextLoc = tempAxTick;

            % Set axis name to turn off and text alignment string
            axName = 'YAxis';
            horizAlignStr = 'right';
            vertAlignStr = 'middle';
    end

    % Plot axis
    hold on
    plot(x,y,'k')
    % Plot tics
    plot(xTick,yTick,'k-')

    % Create text for axis
    if plotTickLabels
        for i = 1:nTick
            text(xTextLoc(i),yTextLoc(i),axisTickLabels{axNo,i}, ...
                'HorizontalAlignment',horizAlignStr, ...
                'VerticalAlignment',vertAlignStr)
        end
    end
end

% Reset axis limits
set(gca,'XLim',xLimNew,'YLim',yLimNew)

% Turn off axis
ax = gca;
ax.(axName).Visible = 'off';
ax.Color = 'none';