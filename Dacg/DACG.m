%-------------------------------------------------------------------------%
%
% File: DACG(x0,f1,K,H,A,L,tol,imax)
%
% Goal: script that performs the DACG method
%
% Inputs: x0:   initial eigenvector estimate 
%               (obtained by some cheap iterations of the power method)
%         f1:   vector containing the function values
%         K:    1/||f1||^2
%         H:    diagonal elements of matrix Theta in LAMBDA x = lambda Theta x
%         A:    preconditioner to accelerate DACG convergence
%         L:    Cholesky Factorization of matrix A
%         tol:  tolerance for the exit test (based on the eigenresidual norm)
%         imax: maximum number of iterations allowed
%
% Calls on: ApplyTheta, ApplyGamma
%
% Outputs: eigenvector: the approximate computed eigenvector
%          ierr:        error flag (set to -1 in case of error)
%                       The test is needed because the matrices involved
%                       could sometimes be (numerically) not strictly 
%                       positive definite. In this case the iteration stops
%                       and the appoximate eigenvector is discarded.
%                       The code recovers by using as approximate eigenvector the
%                       initial approximation provided by the power method. 
%
% Remarks: For more details on the DACG method please Refer to 
%          [L. Bergamaschi, G. Gambolati, G. Pini, 
%          Asymptotic convergence of conjugate gradient methods for the 
%          partial symmetric eigenproblem, Numer. Linear Algebra Appl. 4
%          (1997), pp. 69--84]
%
% Authors: S. De Marchi, A. Martínez, E. Perracchione,
%          Universita' di Padova,
%          Dipartimento di Matematica "Tullio Levi-Civita".
%
% Last modified: 10/10/17.
%
%-------------------------------------------------------------------------%
function [eigenvector,ierr] = DACG(x0,f1,K,H,A,L,tol,imax)
% Initialize
ierr = 0; x = x0;
% Compute first eigenvalue approximation
xB = ApplyTheta(H,x); xA = ApplyGamma(f1,K,L,x); 
gamma = x'*xA; eta = x'*xB;
lambda = gamma/eta;
if lambda<0
    ierr= - 1;
    return;
end
% Compute the first eigenvector approximation and residual vector
eigenvector = x/sqrt(eta);  r = xA - lambda*xB;
res = norm(r)/sqrt(eta);
res0 = res; iter = 1; beta = 0; den0 = 1;
resvec = zeros(imax,1);
p = zeros(length(x),1);
% Start iteration
while (res>tol*res0 & iter<imax)
    % Compute the gradient of the Rayleigh Quotient
    g = r*2/eta; gp = A*g; den = g'*gp; 
    if iter>1
        beta = den/den0;
    end
    gpold = gp; gold = g; den0 = den;
    % Compute new search direction
    p = gp + beta*p; pB = ApplyTheta(H,p);
    pA = ApplyGamma(f1,K,L,p);
    % Compute scalar alpha
    a = p'*xA;b=p'*pA;c = p'*xB;d = p'*pB;
    aa = b*c - a*d; bb = gamma*d - eta*b; cc = eta*a - gamma*c;
    delta = bb^2 - 4*aa*cc; alpha = (bb + sqrt(delta))/(2*aa);
    % Update the eigenvector approximation
    x = x + alpha*p; xA = xA + alpha*pA; xB = xB + alpha*pB;
    % Update the approximate eigenvalue
    gamma = gamma + 2*a*alpha + b*alpha^2;
    eta   = eta   + 2*c*alpha + d*alpha^2;
    lambda = gamma/eta;
    if lambda<0
        ierr = -1;
        return;
    end
    % Compute new residual and residual norm
    r = xA - lambda*xB;
    resvec(iter) = norm(r)/sqrt(eta);
    res = resvec(iter);
    iter = iter + 1;
end
eigenvector = x/sqrt(eta);
iter = iter-1;
