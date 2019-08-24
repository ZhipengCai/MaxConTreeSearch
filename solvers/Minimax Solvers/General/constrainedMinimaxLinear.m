function [sol, fval, idxActive] = constrainedMinimaxLinear(A,b,AEnforce,bEnforce,solInit)

fGurobi = [zeros(size(solInit)); 1];
AGurobi = [[A; -A], -ones(size(A,1)*2,1)];
bGurobi = [b;-b];

AGurobi = [AGurobi ; [[AEnforce; -AEnforce], zeros(2,1)]];
bGurobi = [bGurobi; [bEnforce;-bEnforce]];

[solGurobi, fval] = linprogGurobi(fGurobi, AGurobi, bGurobi);
sol = solGurobi(1:end-1);
res = abs(A*sol-b);

idxActive = find(abs(res-fval)<1e-8);


if numel(idxActive) == 0
    disp(['idxActive is empty, res has' num2str(numel(res)) ' elements; maxres = ' num2str(max(res)) '; fval = ' num2str(fval)]);
    pause;
end
end