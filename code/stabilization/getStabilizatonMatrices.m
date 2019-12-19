function [beta, o] = getStabilizatonMatrices(fa)
% A function to get matrices so latent state (stabilized representation of
% neural signals) can be calculated as:
%
%   l_t = beta*y_t + o,
%
%   where l_t is latent state, y_t is vector of neural data and beta and o
%   are matrices calculated by this function. 
%
% Inputs:
%
%   fa: A FA model with stabilization parameters.  See fitBaseStabilizer
%   for a description of the parameters.
%
% Outputs:
%
%   beta, o: The matrices in the above equation. 
%
% Author: William Bishop, bishopw@janelia.hhmi.org

beta = fa.C'/(fa.C*fa.C' + fa.psi); 
o = -beta*fa.d;


