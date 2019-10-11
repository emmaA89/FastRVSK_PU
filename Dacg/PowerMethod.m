%-------------------------------------------------------------------------%
%
% File: PowerMethod(x0,IM,tol)
%
% Goal: script that performs some iterations of the power method 
%
% Inputs: x0:       the initial vector
%         IM:       the matrix
%         tol:      the tolerance for the exit test
%                   based on the asbolute eigenresidual norm
%
% Outputs: z:       approximate eigenvector corresponding to the 
%                   maximum eigenvalue of IM
%
% Authors: S. De Marchi, A. Martínez, E. Perracchione,
%          Universita' di Padova,
%          Dipartimento di Matematica "Tullio Levi-Civita".
%
% Last modified: 10/10/17.
%
%-------------------------------------------------------------------------%
function z = PowerMethod(x0,IM,tol)
% Initialize
maxiter=10; %Maximum number of iterations allowed
iter = 0;   
normres = 100; 
% Iterative step of the power method
while normres > tol  && iter < maxiter 
    iter = iter + 1;
    z = IM*x0;
    eig =  x0'*z;
    res = (z - eig*x0)/eig;
    normres = norm(res);
    al = norm(z,2);
    z = z/al;
    x0 = z;
end


