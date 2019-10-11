%-------------------------------------------------------------------------%
%
% File: DistanceMatrixCSRBF(dsites,puctrs,wep,M)
%
% Goal: forms the distance matrix of two sets of points for
%       compactly supported RBFs.
%
% Inputs: dsites:   matrix representing a set evaluation points
%         puctrs:   matrix representing a set of PU centres
%         wep:      inverse of CSRBF's support
%         M:        space dimension
%
% Outputs:  DM: NxS SPARSE matrix that contains the Euclidean u-distance
%              (u=1-e*r) between the i-th data site and the j-th center
%              in the i,j position
%
% Calls on: IntegerBased_MD_Structure, IntegerBased_MD_Neighbourhood,
%           IntegerBased_MD_RangeSearch, IntegerBased_MD_ContainingQuery.
%
% Remark:  1. The weight functions are in shifted form
%          2. The original routine is from 
%           [G.E. Fasshauer, Meshfree Approximation Methods with Matlab,
%           World Scientific, Singapore, 2007].
%           Such routine is here modified to work with the integer-based
%           data structure.
%
% Authors: S. De Marchi, A. Martínez, E. Perracchione,
%          Universita' di Padova,
%          Dipartimento di Matematica "Tullio Levi-Civita".
%
% Last modified: 10/10/17.
%
%-------------------------------------------------------------------------%
function DM = DistanceMatrixCSRBF(dsites,puctrs,wep,M)
N = size(dsites,1); S = size(puctrs,1); rbfradius = 1./wep;
nzmax = 25*N; rowidx = zeros(1,nzmax); colidx = zeros(1,nzmax);
validx = zeros(1,nzmax); istart = 1; iend = 0;
% Define the number of square cells in one direction
q = ceil(1./(rbfradius));
 if S > N % Faster if more centres than data sites
    % Build partitioning structure
    idx_ct = IntegerBased_MD_Structure(puctrs,q,rbfradius,M);
    for i = 1:N
        index1 = IntegerBased_MD_ContainingQuery(dsites(i,:),q,...
            rbfradius,M);
        % Find data sites located in the j-th subdomain
        [dxx dx] = IntegerBased_MD_Neighbourhood(puctrs,idx_ct,...
            index1,q,M,1);
        [idx dist] = IntegerBased_MD_RangeSearch(dsites(i,:),...
            rbfradius,dxx,dx);
        newentries = length(idx); iend = iend + newentries;
        rowidx(istart:iend) = repmat(i,1,newentries);
        colidx(istart:iend) = idx'; validx(istart:iend) = 1-wep*dist';
        istart = istart + newentries;
    end
else
    idx_ct = IntegerBased_MD_Structure(dsites,q,rbfradius,M);
     for j = 1:S
         % Find the square cell containing the j-th centre
         index1 = IntegerBased_MD_ContainingQuery(puctrs(j,:),q,...
             rbfradius,M);
         % For each centre, find the data sites in its support
         [dxx dx] = IntegerBased_MD_Neighbourhood(dsites,idx_ct,...
             index1,q,M,1);
         [idx dist] = IntegerBased_MD_RangeSearch(puctrs(j,:),...
             rbfradius,dxx,dx);
         newentries = length(idx); iend = iend + newentries;
         rowidx(istart:iend) = idx';
         colidx(istart:iend) = repmat(j,1,newentries);         
         validx(istart:iend) = 1-wep*dist'; istart = istart + newentries;
     end
 end
idx = find(rowidx); DM = sparse(rowidx(idx),colidx(idx),validx(idx),N,S);