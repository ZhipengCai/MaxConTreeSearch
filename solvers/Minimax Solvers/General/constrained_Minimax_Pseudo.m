function [res, sol, idxActive] = constrained_Minimax_Pseudo(A, b, c, d, AEnforce, bEnforce, cEnforce, dEnforce, th, dim1, solInit, acc)
if nargin <12
    acc = 0;
end
options = optimoptions('fminimax','Display','off', 'ConstraintTolerance', 1e-9, 'FunctionTolerance',  1e-9, 'OptimalityTolerance', 1e-9, 'StepTolerance', 1e-9, 'TolConSQP', 1e-9,'TolConSQP', 1e-9,'MaxIterations',1000); % Minimize abs. values

[sol, res, fval, ~] = fminimax(@resEval, solInit, [], [], [], [], [], [], @mycon, options);

idxActive = find(abs(res-fval)<1e-6);
if numel(idxActive) == 0
    idxActive = find(abs(res-max(res))<1e-6);
end

    function f = resEval(x)
        f = resPseudoConvex(A,b,c,d,x,th,dim1);
    end

    function [c,ceq] = mycon(x)
        c = resPseudoConvex(AEnforce, bEnforce, cEnforce, dEnforce, x, th, dim1)-acc;
        ceq = [];
    end

end