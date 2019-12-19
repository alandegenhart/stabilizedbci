function [postMeans, postCovInds, postCovs] = inferFALatents(faMdl, x, varargin)
% Given an FA model and observed data this function will infer the 
% posterior distribution over latent states. 
%
% Usage: [postMeans, postCovInds, postCovs] = inferFALatents(faMdl, x, varargin)
%
% Inputs: 
%
%   faMdl - A structure with the following fields:
%
%       d - A P by 1 vector of the model mean, where P is the number of
%       observed variables. 
%
%       C - A P by L matrix giving the estimated loading matrix, where L is the
%       number of latent variables. 
%
%       psi - A P by P diagonal matrix of the estimated private noise
%       variances. 
%
%   x - A N by P matrix of data; each row is a sample; columns are variables.
%   Missing values are indicated by NAN. 
%
% Optional Inputs: All optional inputs should be entered in string-value
% pair format. 
%
%   NO_ASSERT - True if assertions should be skipped.  Good for efficiency
%   buy may not catch unexpected inputs.  Default: false. 
%
% Outputs: 
%
%   postMeans - A N by L matrix of the posterior mean of the latent
%   variables for each corresponding sample in x. 
%
%   postCovInds - An array of length N. Each entry is the index into 
%   postCovs (see below) for the posterior covariance matrix for each 
%   corresponding sample in x. 
%
%   postCovs - An L by L by B array, where B is the number of posterior
%   covariance matrices, containing the posterior covariance matrices for
%   the posterior distribution over latents. 
%
% Author: William Bishop, bishopw@janelia.hhmi.org

NO_ASSERT = false;
warnOpts(assignOpts(varargin)); 

if ~NO_ASSERT
    assertFAMdl(faMdl); 
end

nSmps = size(x,1); 
nLatents = size(faMdl.C,2); 

% Subtract the mean from x
xCtr = bsxfun(@minus, x, faMdl.d');

% Arrange centered data into blocks according to what variables are observed
xCtrBlocked = smpMat2BlockStr(xCtr);

% Estimate latents
latentStr = fastFAEStep(xCtrBlocked, faMdl.C, faMdl.psi);

% Put latents back into a big matrix 
nBlks = length(latentStr); 
postMeans = nan(nSmps,nLatents); 
postCovInds = nan(nSmps,1); 
postCovs = nan(nLatents, nLatents, nBlks); 

for bI = 1:nBlks
    postMeans(xCtrBlocked(bI).origRows,:) = latentStr(bI).postMeans;
    postCovInds(xCtrBlocked(bI).origRows) = bI;
    postCovs(:,:,bI) = latentStr(bI).postCov; 
end