function warnOpts(optsCellIn)
% This is a function that is intended to work in conjuction with
% assignOpts.  It will accept a cell of string-value option arguments in
% the same form as output by that function (see assignOpts for more
% details) and for each string it will print a warning message indicating
% that a variable with that name was not assigned a value. 
%
% This functionis is based on one in code provided to Will Bishop by the 
% Shenoy lab at Stanford University during his time at JHU/APL. 
%
% Usage: warnOpts(optsCellIn)
%
% Inputs: optsCellIn - a cell of string-value pair option arguments to
%   print warning messages for. 
%
% Outputs: None. 
%
% Author: William Bishop, bishopw@janelia.hhmi.org

lOptsCell = length(optsCellIn); 

% Make sure the length of optsCell is even
assert(mod(lOptsCell,2) == 0, 'Option names/values must come in pairs.'); 

% See how many options we have to assign
nOpts = lOptsCell/2; 

% For each option, print a warning
for o = 1:nOpts
    curVarInd = (o-1)*2+1;
    curVarName = optsCellIn{curVarInd};    
    warning('CUSTOM_UTILS:unassigned', ['The variable ', curVarName, ' was not assigned.']); 
end