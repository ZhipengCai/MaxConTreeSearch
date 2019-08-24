function [y] = solveY_BCO(s,delta)
    [sSorted, idxSorted] = sort(s,'ascend');
    N = length(s);
    y = [ones(N,1)];
    y(idxSorted(delta+find(sSorted(delta+1:N)>0)))=0;
end