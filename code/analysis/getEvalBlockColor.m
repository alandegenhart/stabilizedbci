function [cPatch,cSolid] = getEvalBlockColor(blockName)
% getEvalBlockColor         Get color for evaluation block

switch lower(blockName)
    case 'baseline evaluation'
        % Blue
        cPatch = [210 230 255]/255;
        cSolid = [0 0 1];
    case 'stitching evaluation'
        % Green
        cPatch = [210 255 210]/255;
        cSolid = [14 135 0]/255;
    case 'perturbation evaluation'
        % Red
        cPatch = [255 220 210]/255;
        cSolid = [1 0 0];
    case {'baseline washout','baseline'}
        % Gray
        cPatch = ones(1,3) * 0.85;
        cSolid = [0 0 0];
end
    