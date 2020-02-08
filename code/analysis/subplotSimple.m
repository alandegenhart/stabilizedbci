% subplotSimple  Generate subplot with specific sizing.
%
% This function replicates the functionality of the 'subplot' command, but
% allows for additional functionality. In particular, the optional argument
% 'Ax' can be used in order to automatically determine the subplot location
% based on those retured by the 'calcFigureSize' function.
%
% Usage:
%   [h] = subplotSimple(nRow, nCol, plotNo)
%
% Inputs:
%   nRow    Number of rows
%   nCol    Number of columns
%   plotNo  Current plot number
%
% Optional inputs:
%   Ax      Axis structure returned by calcFigureSize
%   xMarg   Array containing left and right margins (ignored if Ax is
%           specified)
%   yMarg   Array containing bottom and top margins (ignored if Ax is
%           specified)
%   xSpace  Inter-plot x-spacing (ignored if Ax is specified)
%   ySpace  Inter-plot y-spacing (ignored if Ax is specified)
%
% Outputs:
%    h      Axis handle
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

function [h] = subplotSimple(nRow, nCol, plotNo, varargin)

% Optional arguments
Ax = struct();
xMarg = [.1 .1];  % L, R margin
yMarg = [.1 .1];   % B, T margin
xSpace = 0;         % Inter-plot spacing (x)
ySpace = 0;         % Inter-plot spacing (y)

% Assign optional arguments
assignOpts(varargin);

% If 'Ax' structure is provided as an input, set values accordingly
if ~isempty(Ax)
    xMarg = Ax.xMarg;
    yMarg = Ax.yMarg;
    xSpace = Ax.xSp;
    ySpace = Ax.ySp;
end

if nRow == 1
    ySpace = 0;
end

if nCol == 1
    xSpace = 0;
end

w = (1 - sum(xMarg) - (nCol-1)*xSpace) / nCol;
h = (1 - sum(yMarg) - (nRow-1)*ySpace) / nRow;

% Find positions for each plot
xPos = xMarg(1):(w+xSpace):(1-xMarg(2));
yPos = (1-yMarg(2)-h):-(h+ySpace):yMarg(1); 
xPos = xPos(1:nCol);

% Handle case where only one row is to be plotted
if nRow == 1
    yPos = yMarg(1);
else
    yPos = yPos(1:nRow);
end

% Convert plot index into row/column index
rowIdx = ceil(plotNo/nCol);
colIdx = rem(plotNo,nCol);
colIdx(colIdx == 0) = nCol; % Handle case for the last column (remainder is 0)

% Handle case for if there are multiple inputs for the plotNo
% Determine the starting position
xPosStart = min(xPos(colIdx));
yPosStart = min(yPos(rowIdx));

% Determine the end edge
xPosEdge = max(xPos(colIdx))+w;
yPosEdge = max(yPos(rowIdx))+h;

% Redefine w and h if necessary
w = xPosEdge - xPosStart;
h = yPosEdge - yPosStart;

% Creat plot
h = subplot('Position',[xPosStart,yPosStart,w,h]);
