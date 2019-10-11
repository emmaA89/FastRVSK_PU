function gridpoints = MakeSDGrid(s,neval)

% The function is from
% [G.E. Fasshauer, Meshfree Approximation Methods with Matlab,
% World Scientific, Singapore, 2007]

if (s==1)
    gridpoints = linspace(0,1,neval)';
    return;
end
outputarg = 'x1';
for d = 2:s
    outputarg = strcat(outputarg,',x',int2str(d));
end
makegrid = strcat('[',outputarg,'] = ndgrid(linspace(0,1,neval));');
eval(makegrid);
gridpoints = zeros(neval^s,s);
for d = 1:s
    matrices = strcat('gridpoints(:,d) = x',int2str(d),'(:);');
    eval(matrices);
end
