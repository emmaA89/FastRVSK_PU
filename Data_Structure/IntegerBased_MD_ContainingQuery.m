function [index] = IntegerBased_MD_ContainingQuery(puctr,q,puradius,M)

% The function is from
% http://hdl.handle.net/2318/1559094
%
% Remarks: Refer also to 
%          [R. Cavoretto, A. De Rossi, E. Perracchione,
%          Optimal selection of local approximants in RBF-PU interpolation, 
%          to appear on J. Sci. Comput. (2017)]

idx = ceil(puctr./puradius); k = 1:M-1;
idx(idx == 0) = 1;
index = sum((idx(k)-1).*q.^(M-k)) + idx(end);