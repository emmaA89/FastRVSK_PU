%-------------------------------------------------------------------------%
%
% File: DistanceMatrix_VSK(dsites,ctrs,h,puctrs,puradius)
%
% Goal: forms the distance matrix of two sets of points via VSKs
%
% Inputs: dsites:   matrix representing a set data points
%         ctrs:     matrix representing a set of RBF centres
%         h:        parameter for the scale function
%         puctrs:   the PU centre
%
% Outputs:  DM: matrix that contains the Euclidean distances via VSKs
%
% Calls on: Scale_Function
%
% Remarks: Refer also to 
%          [M. Bozzini, L. Lenarduzzi, M. Rossini, R. Schaback, 
%          Interpolation with variably scaled kernels, 
%          IMA J. Numer. Anal. 35 (2015), 199--219]
%
% Authors: S. De Marchi, A. Martínez, E. Perracchione,
%          Universita' di Padova,
%          Dipartimento di Matematica "Tullio Levi-Civita".
%
% Last modified: 10/10/17.
%
%-------------------------------------------------------------------------%
function DM = DistanceMatrix_VSK(dsites,ctrs,h,puctrs)
[M,s] = size(dsites); [N,s] = size(ctrs); % Initialize
DM = zeros(M,N); s = s+1; % Increment the dimension
dsites = [dsites Scale_Function(dsites,h,puctrs)'];
ctrs = [ctrs Scale_Function(ctrs,h,puctrs)'];
for d=1:s
    [dr,cc] = ndgrid(dsites(:,d),ctrs(:,d));
    DM = DM + (dr-cc).^2;
end
DM = sqrt(DM);