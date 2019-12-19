function T = learnOptimalOrthonormalTransformation(m1, m2)
% This is function to learn the optimal orthogonal transformation between 
% two matrices such that 
%
%   ||m1 - m2*T'||_F 
%
% is minimized, where m1 and m2 are two matrices and T is the learned
% transformation and T' indicates the transpose of T. 
%
% This is based on the paper "A Generalized Solution to the Orthogonal
% Procrustes Problem" by Peter H. Schonemann. 
%
% Usage: T = learnOptimalOrthonormalTransformation(m1, m2)
%
% Input: 
%
%   m1, m2 - the two matrices to learn the transformation for. 
%
% Output:
%
%   T - the learned transformation. 
%
% Author: William Bishop, bishopw@janelia.hhmi.org

S = m1'*m2; 
[U, ~, V] = svd(S); 
T = U*V'; 
