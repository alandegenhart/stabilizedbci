function [fW,fH,Ax] = calcFigureSize(nRow,nCol,axW,axH,axSp,varargin)
%
% Inputs:
%   nRow        Number of rows
%   nCol        Number of columns
%   axW         Axis width (pixels)
%   axH         Axis height (pixels)
%   sp          Axis spacing (x and y) (pixels)

% Optional arguments

margin = 75; % Margin (px)

assignOpts(varargin);

% Calculate figure size
fW = axW*nCol + axSp*(nCol-1) + margin*2;
fH = axH*nRow + axSp*(nRow-1) + margin*2;

% Calculate fractional values (for use with the 'subplotSimple' command
Ax.xMarg = margin/fW;
Ax.xSp = axSp/fW;
Ax.yMarg = margin/fH;
Ax.ySp = axSp/fH;
Ax.nRow = nRow;
Ax.nCol = nCol;