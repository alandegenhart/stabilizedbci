% getEvalBlockColor  Get color for specified evaluation block.
%
% This function returns RGB color codes for the provided block name.
%
% Usage:
%   [cPatch, cSolid] = getEvalBlockColor(blockName)
%
% Inputs:
%   blockName   Name of block to get color information for
%
% Outputs:
%   cPatch      Light color (used for background/fill)
%   cSolid      Dark color (used for solid lines)
%
% @ Alan Degenhart -- alan.degenhart@gmail.com

function [cPatch,cSolid] = getEvalBlockColor(blockName)

switch lower(blockName)
    case 'baseline evaluation'
        % Blue
        cPatch = [210, 230, 255]/255;
        cSolid = [0, 0, 1];
    case 'stitching evaluation'
        % Green
        cPatch = [210, 255, 210]/255;
        cSolid = [14, 135, 0]/255;
    case 'perturbation evaluation'
        % Red
        cPatch = [255, 220, 210]/255;
        cSolid = [1, 0, 0];
    case {'baseline washout','baseline'}
        % Purple
        cPatch = [250, 155, 236]/255;
        cSolid = [181, 0, 154]/255;
    case 'stitching'
        % Gray
        cPatch = ones(1,3) * 0.9;
        cSolid = [0, 0, 0];
    otherwise
        error('Invalid block name specified.')
end
    