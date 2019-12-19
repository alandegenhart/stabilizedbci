function assertFAMdl(mld)
% This is a function to assert that an FA model is valid. 
%
% Usage: assertValidFAMdl(mld)
%
% Inputs: 
%
%   mdl - a structure of the FA model, with fields C, d, psi.  See the
%   function fitFAWithMissingVls.m for more information. 
%
% Outputs: 
%
%   None. 
%
% Author: William Bishop, bishopw@janelia.hhmi.org

% Make sure all parameters are real
assert(all(all(isreal(mld.C))), 'C must contain only real numbers.'); 
assert(all(isreal(mld.d)), 'd must contain only real numbers.');
assert(all(all(isreal(mld.psi))), 'psi must contain only real numbers.');

% Check dimensions
nDObs = length(mld.d); 
assert(length(mld.psi) == nDObs, 'Number of observations between psi and d is inconsistent.');
if size(mld.C,2) >= 0
    assert(size(mld.C,1) == nDObs, 'Number of observations between C and d is inconsistent.'); 
end

% Make sure psi is a legitimate diagonal covariance matrix 
diagPsi = diag(mld.psi);
assert(ismatrix(mld.psi), 'psi must be a matrix.'); 
assert(all(all(mld.psi == diag(diagPsi))), 'psi must be diagonal'); 
assert(all(diagPsi >= 0), 'All entries in psi must be greater than or equal to 0.'); 

