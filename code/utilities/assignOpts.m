function optsCellOut = assignOpts(optsCellIn)
% This is a function that is based on code provided to Will Bishop
% by the Shenoy lab at Stanford University during his time at JHU/APL. 
%
% It will accept a cell as input. It is expected the length of this cell 
% will be even and will consist of pairs of entries.  The first entry in 
% each pair should be a variable name and the second entry should be the 
% value for that variable.  
%
% When this function is called it will look through the variables in the
% callers workspace, and if one of the variables in the callers workspace
% has the same name as one of the variables in the input cell, it will
% assign the value in the input cell for that variable to the variable in
% the caller's workspace.  Note that it will also do a check to make sure
% that the value in the input cell is of the same class as the existing
% variable in the caller's workspace. 
%
% If there are any variables in the input cell that are not found in the
% caller's workspace, these will be passed out of the function along with
% their values in the output cell. 
%
% Usage: optsCellOut = assignOpts(optsCellIn)
%
% Inputs: 
%
%   optsCellIn - a cell of variable name - variable value pairs.  The
%   variable names should be candidate variables to look for in the
%   caller's workspace and assign new values to.  For example, the cell:
%
%       {'var1', 1, 'var2', 'cat'} 
%
%   would cause this function to look for the variable named var1 in the
%   caller's workspace and assign it value 1 if it was found.  Similiary,
%   the variable 'var2' would be assigned 'cat' if it was found.  
%
% Outputs: 
%
%   optsCellOut - a cell of variable name - variable value pairs that were
%   not found in the caller's workspace. 
%
% Author: William Bishop, bishopw@janelia.hhmi.org

lOptsCell = length(optsCellIn); 

% Make sure the length of optsCell is even
assert(mod(lOptsCell,2) == 0, 'Option names/values must come in pairs.'); 

% See how many options we have to assign
nOpts = lOptsCell/2; 

% Get a list of variables in the caller's workspace
callerVars = evalin('caller', 'whos');
callerVarNames = {callerVars.name};

% Try to assign each option, keeping track of those which were not assigned
optsNotAssigned = false(1, lOptsCell); 
for o = 1:nOpts
    curVarInd = (o-1)*2+1;
    curValInd = (o-1)*2+2; 
    
    curVarName = optsCellIn{curVarInd};
    curVarVal = optsCellIn{curValInd};
    
    if any(strcmpi(curVarName, callerVarNames))
        
        % Make sure we are not switching classes
        curValClass = class(curVarVal);
        existingValClass = evalin('caller', ['class(', curVarName, ')']); 
        if ~(strcmp(curValClass, 'matlab.graphics.axis.Axes') || strcmp(curValClass, 'matlab.ui.container.Panel'))
            assert(strcmpi(curValClass, existingValClass), ['Variable ', curVarName, ' must be of type ', existingValClass, '.']); 
        end
        % Assuming we passed the last assert, assign the new value
        assignin('caller', curVarName, curVarVal); 
        
    else
        optsNotAssigned([curVarInd, curValInd]) = true; 
    end
end

% Place any unassigned options and their values in the optsCellOut cell
optsCellOut = optsCellIn(optsNotAssigned); 