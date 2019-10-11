function [idx, dist] = IntegerBased_MD_RangeSearch(puctr,puradius,...
    dsites,index)

% The function is from
% http://hdl.handle.net/2318/1559094
%
% Remarks: Refer also to 
%          [R. Cavoretto, A. De Rossi, E. Perracchione,
%          Optimal selection of local approximants in RBF-PU interpolation, 
%          to appear on J. Sci. Comput. (2017)]

N = size(dsites,1); dist = []; idx = []; % Initialize
for i = 1:N
    dist1(i) = norm(puctr - dsites(i,:));
end
if N > 0
    [sort_dist,IX] = sort(dist1);
    N1 = size(sort_dist,2); j1 = 1; j2 = 1; 
    if nargin == 3
        idx = IX; dist = dist1;
    else
        while (j2 <= N1) && (sort_dist(j2) <= puradius)
            idx(j1) = index(IX(j2)); dist(j1) = dist1(IX(j2));
            j1 = j1 + 1; j2 = j2 + 1;
        end
    end
end