%-------------------------------------------------------------------------%
%
% File: ApplyTheta(H,x)
%
% Goal: script that computes the product of matrix Theta of the 
%       generalized eigenvalue problem LAMBDA x = lambda Thetha x
%
% Inputs: H:     vector containing the diagonal elements of matrix Theta
%         x:     vector 
%
% Outputs: ax: the product of Theta times vector x
%
% Authors: S. De Marchi, A. Martínez, E. Perracchione,
%          Universita' di Padova,
%          Dipartimento di Matematica "Tullio Levi-Civita".
%
% Last modified: 10/10/17.
%
%-------------------------------------------------------------------------%
function ax = ApplyTheta(H,x)
ax = H.*x;



