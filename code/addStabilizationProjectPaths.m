function addStabilizationProjectPaths(baseDir)
% Adds all required paths for running the stabilization project code.
%
% Author: William Bishop, bishopw@janelia.hhmi.org

if nargin < 1
    baseDir = pwd;
end

addpath(genpath(baseDir));
