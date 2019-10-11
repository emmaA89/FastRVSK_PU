%-------------------------------------------------------------------------%
%
% File: PU_RBF_VSK_RATIONAL_DAGC(M,dsites,neval,npu,rbf,wf,f,h)
%
% Goal: script that performs partition of unity rational interpolation
%       via VSKs and DACG
%
% Inputs:  M:            space dimension
%          dsites:       NXM matrix representing a set of N data sites
%          neval:        number of evaluation points in one direction
%          npu:          number of PU subdomains in one direction
%          rbf:          radial basis function
%          wf:           weight function
%          f:            the test function
%          h:            parameter for the scale function
%
% Outputs:  epoints: the evaluation points
%           Pf:      the interpolant computed at the evaluation points
%
% Calls on: IntegerBased_MD_Structure, IntegerBased_MD_Neighbourhood,
%           IntegerBased_MD_RangeSearch, IntegerBased_MD_ContainingQuery,
%           DistanceMatrixCSRBF, MakeSDGrid, DistanceMatrix_VSK, 
%           PowerMethod and DACG.
%
% References: 1. [L. Bergamaschi, G. Gambolati, G. Pini, 
%          Asymptotic convergence of conjugate gradient methods for the 
%          partial symmetric eigenproblem, Numer. Linear Algebra Appl. 4
%          (1997), 69--84]
%             2. [M. Bozzini, L. Lenarduzzi, M. Rossini, R. Schaback, 
%          Interpolation with variably scaled kernels, 
%          IMA J. Numer. Anal. 35 (2015), 199--219]
%             3. [R. Cavoretto, A. De Rossi, E. Perracchione,
%          Optimal selection of local approximants in RBF-PU interpolation, 
%          to appear on J. Sci. Comput. (2017)]
%             4. [G. E. Fasshauer, Meshfree approximation methods with 
%          Matlab, World Scientific, Singapore, 2007].
%             5. [S. Jakobsson, B. Andersson, F. Edelvik, 
%          Rational radial basis function interpolation with applications 
%          to antenna design}, J. Comput. Appl. Math. 233 (2009), 889--904.
%
% Authors: S. De Marchi, A. Martínez, E. Perracchione,
%          Universita' di Padova,
%          Dipartimento di Matematica "Tullio Levi-Civita".
%
% Last modified: 10/10/17.
%
%-------------------------------------------------------------------------%
function [epoints Pf] = PU_RBF_VSK_RATIONAL_DAGC(M,dsites,neval,npu,rbf,...
    wf,f,h)
% Create npu^M equally spaced PU centres
puctrs = MakeSDGrid(M,npu);
% Create neval^M equally spaced evaluation points
epoints = MakeSDGrid(M,neval);
puradius = 1./npu;  % Define the PU radius
wep = 1./puradius;  % Parameter for weight function
npu_M = size(puctrs,1); neval_M = size(epoints,1); % Initialize;
rbfctrs = dsites; % Define the RBF centres
Pf = zeros(neval_M,1);  % Initialize
% Compute Shepard evaluation matrix
DM_eval = DistanceMatrixCSRBF(epoints,puctrs,wep,M);
SEM = wf(wep,DM_eval);
SEM = spdiags(1./(SEM*ones(npu_M,1)),0,neval_M,neval_M)*SEM;
% Parameter for integer-based partitioning structure
q = ceil(1./puradius);
% Build the partitioning structure for data sites and evaluation points
idx_ds = IntegerBased_MD_Structure(dsites,q,puradius,M);
idx_ep = IntegerBased_MD_Structure(epoints,q,puradius,M);
CN = []; % Initialize
for j = 1:npu_M
    % Find the block containing the j-th subdomain centre
    index1 = IntegerBased_MD_ContainingQuery(puctrs(j,:),q,puradius,M);
    % Find data sites located in the j-th subdomain
    [dxx dx] = IntegerBased_MD_Neighbourhood(dsites,idx_ds,index1,q,M,1);
    idx = IntegerBased_MD_RangeSearch(puctrs(j,:),puradius,dxx,dx);
    % Build local VSK interpolation matrix for the j-th subdomain
    DM_data = DistanceMatrix_VSK(dsites(idx,:),rbfctrs(idx,:),h,...
        puctrs(j,:));
    IM = rbf(1,DM_data);
    [edxx edx] = IntegerBased_MD_Neighbourhood(epoints,idx_ep,index1,q,...
        M,1);
    CN(j) = condest(IM); % Compute the condition number
    % Perform rational VSK interpolation
    Nr = size(dsites(idx,:),1); f1 = f(dsites(idx,:)); % Initialize
    mu = 10^(-14); B = IM+mu.*eye(Nr); % Diagonal correction
    L = chol(B,'lower'); % Factorize
    % Power method to approximate x0
    x0 = ones(size(IM,1),1); % Initialize
    x0 = PowerMethod(x0,IM,1e-3);
    % perform DACG method
    K = 1.0/sum(f1.^2); H = K*f1.^2 + 1; % Define the matrices
    prec = IM; % Preconditioning
    toldacg = 1e-2; % Tolerance for DACG
    imax=30; % Maximum number of iterations
    [q2,ierr] = DACG(x0,f1,K,H,prec,L,toldacg,imax);
    if ierr==-1
        q2=x0;
    end
    % Solve the linear systems
    T0 = [f1.*q2 q2];
    T1 = L'\(L\T0);
    pAlpha = T1(:,1);
    qAlpha = T1(:,2);
    % Find the evaluation points located in the j-th subdomain
    eidx = IntegerBased_MD_RangeSearch(puctrs(j,:),puradius,edxx,edx);
    if  (~isempty(eidx))
        % Compute local VSK evaluation matrix and the local RBF interpolant
        DM_eval = DistanceMatrix_VSK(epoints(eidx,:),rbfctrs(idx,:),...
            h,puctrs(j,:));
        EM = rbf(1,DM_eval); localfit = (EM*pAlpha)./(EM*qAlpha);
        % Accumulate global fit
        Pf(eidx) = Pf(eidx) + localfit.*SEM(eidx,j);
    end
end
% Compute exact solution
exact = f(epoints);
% Compute errors on evaluation grid and the maximum condition numbers
rms_err = norm(Pf - exact)/sqrt(neval_M);
MCN = max(CN(:));
fprintf('RMS error:         %e\n', rms_err);
fprintf('average_Cond:      %e\n', MCN);