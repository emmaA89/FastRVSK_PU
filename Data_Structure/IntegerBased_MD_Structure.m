function [idx_dsites_k] = IntegerBased_MD_Structure(dsites,q,puradius,M)

% The function is from
% http://hdl.handle.net/2318/1559094
%
% Remarks: Refer also to 
%          [R. Cavoretto, A. De Rossi, E. Perracchione,
%          Optimal selection of local approximants in RBF-PU interpolation, 
%          to appear on J. Sci. Comput. (2017)]


N = size(dsites,1); idx_dsites_k = cell(q^(2*M),1); k = 1:M-1;
for i = 1:N
    idx = ceil(dsites(i,:)./puradius);
    idx(idx == 0) = 1;
    index = sum((idx(k)-1).*q.^(M-k)) + idx(end);
    idx_dsites_k{index} = [idx_dsites_k{index}; i];
end