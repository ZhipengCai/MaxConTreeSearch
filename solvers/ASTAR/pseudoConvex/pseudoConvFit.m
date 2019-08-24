function [sol, outl, UniqueNodeNumber, hInit, ubInit, levelMax, NODIBP, runtime] = pseudoConvFit(A,b,c,d,solInit,th,method,runtimeLimit)
if nargin < 6
    runtimeLimit = -1;
end
tic;
if strcmp(method, 'A*')
    [sol, outl, ~,UniqueNodeNumber,hInit, ubInit, levelMax] = ASTARWithH2_pseudo(A,b,c,d, solInit, th, 1, 0, runtimeLimit,0);
    NODIBP = 0;
elseif strcmp(method, 'A*-NAPA')
    [sol, outl, ~,UniqueNodeNumber,hInit, ubInit, levelMax] = ASTARWithH2_pseudo(A,b,c,d, solInit,  th, 1, 0,runtimeLimit);
    NODIBP = 0;
elseif strcmp(method, 'A*-TOD')
    [sol, outl, ~,UniqueNodeNumber,hInit, ubInit, levelMax, ~, NODIBP] = ASTARWithH2_GORE_pseudo(A,b,c,d, solInit, th, 1, 0, runtimeLimit,0);
elseif strcmp(method, 'A*-NAPA-TOD')    
    [sol, outl, ~,UniqueNodeNumber,hInit, ubInit, levelMax, ~, NODIBP] = ASTARWithH2_GORE_pseudo(A,b,c,d, solInit, th, 1, 0, runtimeLimit);
elseif strcmp(method, 'A*-NAPA-DIBP')    
     [sol, outl, ~,UniqueNodeNumber,hInit, ubInit, levelMax, ~, NODIBP] = ASTARWithH2_DIBP_pseudoV2(A,b,c,d, solInit, th, 1, 0, runtimeLimit);
end
runtime = toc;
end
