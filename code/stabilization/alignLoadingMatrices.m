function [W, alignChs] = alignLoadingMatrices(m1, m2, n, th)
% This is a function to align the loading matrices of two FA models.
%
% Usage: W = alignLoadingMatrices(m1, m2, n, th)
%
% Inputs:
%
%   m1: the loading matrix for the original model
%
%   m2: the loading matrix for the new model
%
%   n: the number of rows to use for alignment
%
%   th: the threshold to use when screening out rows for possibly
%   alignment.  Any row which has an l_2 norm less than th in either m_1 or
%   m_2 will not be considered when selecting rows to use for alignment.
%
% Outputs:
%
%   W: the matrix that will align m2 to m1 so that m1 - m2*W' is minimized.
%
%   alignChs: rows/cols of m1 and m2 that were used for alignment
%
%
% Author: William Bishop, bishopw@janelia.hhmi.org

alignChs = identifyStableLoadingRows(m1, m2, n, th);
W = learnOptimalOrthonormalTransformation(m1(alignChs,:), m2(alignChs,:));  



