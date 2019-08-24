function [res, sol, idxActive] = constrainedMinimaxLinearV2(A, b, ACon, bCon, solInit)

fGurobi = [zeros(size(solInit)); 1];
AGurobi = [[A; -A], -ones(size(A,1)*2, 1)];
bGurobi = [b;-b];

AGurobi = [AGurobi ; [ACon, zeros(size(ACon,1),1)]];
bGurobi = [bGurobi; bCon];

options = optimoptions('linprog', 'Display', 'none');

[solGurobi, fval, exitflag] = linprog(fGurobi, AGurobi, bGurobi, [], [], [], [], options);

if exitflag == -2
    sol = solInit;
    idxActive = [];
    res = [];
    return;
end

sol = solGurobi(1:end-1);
res = abs(A*sol-b);

idxActive = find(abs(res-fval)<1e-8);

% % used to test whether the active set is empty
% if numel(idxActive) == 0
%     disp(['idxActive is empty, res has' num2str(numel(res)) ' elements; maxres = ' num2str(max(res)) '; fval = ' num2str(fval)]);
%     pause;
% end

end