function [dxx dx] = IntegerBased_MD_Neighbourhood(dsites,idx_ds,index1,...
    q,M,t)

% The function is from
% http://hdl.handle.net/2318/1559094
%
% Remarks: Refer also to 
%          [R. Cavoretto, A. De Rossi, E. Perracchione,
%          Optimal selection of local approximants in RBF-PU interpolation, 
%          to appear on J. Sci. Comput. (2017)]

neigh = []; k = M-1; index = index1; k1 = 1:t; % Initialize
while k  > 0
     neigh = [index1+k1.*q.^k,index1-k1.*q.^k];
    if k - 1 > 0
        neigh = [neigh,neigh+q.^(k-1),neigh-q.^(k-1)];
    end
    k = k - 1;
end
k2 = 1:t; k3 = k2;
for i = 1:length(neigh)
    neighplus(k2) = neigh(i) + k3;
    neighminus(k2) = neigh(i) - k3;
    k2 = k2 + t;
end
neigh = [neigh,index1+k1,index1-k1,neighplus,neighminus];
j = find(neigh > 0 & neigh <= q^M);
index = [index; neigh(j)'];
dxx = []; dx = []; 
for p = 1:length(index)
    dxx = [dxx;dsites(idx_ds{index(p)},:)];
    dx = [dx;idx_ds{index(p)}];
end