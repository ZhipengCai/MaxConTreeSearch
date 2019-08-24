function [res, sol, idxActive] = constrained_Minimax_Linear(A, b, AEnforce, bEnforce, th, solInit, fval)
if nargin < 7
    fval = 0;
end
options = optimoptions('linprog', 'Display', 'none');

fGurobi = [zeros(size(solInit)); 1];

AGurobi = [[A; -A], -ones(size(A,1)*2,1)];
bGurobi = [b;-b];

AGurobi = [AGurobi ; [[AEnforce; -AEnforce], zeros(size(AEnforce,1)*2,1)]];
bGurobi = [bGurobi; [bEnforce;-bEnforce]+th];

[solGurobi, fval, exitflag] = linprog(fGurobi, AGurobi, bGurobi, [], [], [], [], options);

sol = solGurobi(1:end-1);
res = abs(A*sol-b);

idxActive = find(abs(res-abs(fval))<1e-8);

if numel(idxActive) == 0
    idxActive = find(abs(res-max(res))<1e-8);
    disp(['idxActive is empty, res has' num2str(numel(res)) ' elements; maxres = ' num2str(max(res)) '; fval = ' num2str(fval), '; exitflag = ' numel(exitflag)]);
end

end