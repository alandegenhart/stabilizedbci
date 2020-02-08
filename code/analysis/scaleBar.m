% scaleBar  Add scale bar with units to plot.
%
% Usage:
% scaleBar(ax, d, uniStr)
%
% Inputs:
%   ax      Plot axis handle
%   d       Length of scale bar
%   uniStr  Unit string (e.g., 'mm')
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

function scaleBar(ax, d, unitStr, varargin)

% Handle optional arguments
xScale = true;
yScale = false;
fontSize = 10;

assignOpts(varargin);

% Get axis limits (used for placing scale bar)
xLim = get(ax,'XLim');
yLim = get(ax,'YLim');

% Get lower-right point for scale bar.  This function could be updated in
% the future to make the location an optional input argument.
x = xLim(1) + (xLim(2) - xLim(1))*.95;
y = yLim(1) + (yLim(2) - yLim(1))*.05;

% Plot x scale bar
if xScale
    plot([x-d x],[y y],'k-');
    
    % Place scale text
    h = text(x-d/2,y,sprintf('%d %s',d,unitStr), 'FontSize', fontSize);
    set(h,'HorizontalAlignment','center','VerticalAlignment','bottom')
end

% Plot y scale bar
if yScale
    plot([x x],[y y+d],'k-'); %#ok<*UNRCH>
    
    % Place scale text
    h = text(x,y+d/2,sprintf('%d %s',d,unitStr), 'FontSize', fontSize);
    set(h,'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'rotation',90)
end