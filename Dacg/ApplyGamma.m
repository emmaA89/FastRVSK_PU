%-------------------------------------------------------------------------%
%
% File: ApplyGamma(f1,K,L,x) 
%
% Goal: script that computes the product of matrix LAMBDA of the
%       generalized eigenvalue problem LAMBDA x = lambda Theta x, 
%       times vector x
%       LAMBDA = 1/||f1||^2 D^T A^-1 D + A^-1
%       D: diagonal matrix with the function values as diagonal elements
%
% Inputs: f1:    vector containing the function values
%         K:     1/||f1||^2
%         L:     Cholesky factorization of each local interpolation matrix A (IM in the code)
%         x:     vector
%
% Outputs: ax: the product of matrix LAMBDA times vector x,
%              (1/||f1||^2 D^T A^-1 D + A^-1)x
%
% Authors: S. De Marchi, A. Martínez, E. Perracchione,
%          Universita' di Padova,
%          Dipartimento di Matematica "Tullio Levi-Civita".
%
% Last modified: 10/10/17.
%
%-------------------------------------------------------------------------%
function ax = ApplyGamma(f1,K,L,x) 
v = f1.*x; T0 = [x v]; T1 = L'\(L\T0);
ax = K*f1.*T1(:,2)+T1(:,1);
