function AX = subplotSimple(nRow,nCol,plotNo,varargin)
% [AX] = subplotSimple(nRow,nCol,plotNo)
%
% Generate subplot with no whitespace.
%
% This function replicates the functionality of the 'subplot' command, but
% generates plots with no whitespace (i.e., there is no distance between
% plots).
%
% Inputs:
%   nRow    Number of rows
%   nCol    Number of columns
%
% Outputs:
%   AX      Axes handle
%
% Author:       Alan D. Degenhart
% Date Created: 2016.05.21
% Last Updated: 2016.05.21

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
    xMarg = [Ax.xMarg Ax.xMarg];
    yMarg = [Ax.yMarg Ax.yMarg];
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
AX = subplot('Position',[xPosStart,yPosStart,w,h]);

% if colIdx == 0 % Handle case for the last column (remainder is 0)
%     colIdx = nCol;
% end
% % Creat plot
% AX = subplot('Position',[xPos(colIdx),yPos(rowIdx),w,h]);