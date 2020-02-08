% calcFigureSize  Get figure and subplot data for specified subplot
% configuration.
%
% This function determines the figure and subplot size information for the
% desired plot configuration.
%
% Usage
%   [fW, fH, Ax] = calcFigureSize(nRow, nCol, axW, axH, axSp)
%
% Inputs:
%   nRow        Number of rows
%   nCol        Number of columns
%   axW         Axis width (pixels)
%   axH         Axis height (pixels)
%   axSp        Axis spacing (x and y) (pixels)
%
% Optional arguments:
%   xMargin     Margins along x-axis of figure (pixels)
%   yMargin     Margins along y-axis of figure (pixels)
%
% Outputs:
%   fW          Figure width (pixels)
%   fH          Figure height (pixels)
%   Ax          Structure containing information for creating subplots.
%               Used along with the subplotSimple command.
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

function [fW, fH, Ax] = calcFigureSize(nRow, nCol, axW, axH, axSp, varargin)


% Optional arguments
xMargin = [75, 75];
yMargin = [75, 75]; % Margin (px)

assignOpts(varargin);

% Calculate figure size
fW = axW*nCol + axSp*(nCol-1) + sum(xMargin);
fH = axH*nRow + axSp*(nRow-1) + sum(yMargin);

% Calculate fractional values (for use with the 'subplotSimple' command
Ax.fW = fW;
Ax.fH = fH;
Ax.axW = axW;
Ax.axH = axH;
Ax.axSp = axSp;
Ax.xMarg = xMargin/fW;
Ax.xSp = axSp/fW;
Ax.yMarg = yMargin/fH;
Ax.ySp = axSp/fH;
Ax.nRow = nRow;
Ax.nCol = nCol;