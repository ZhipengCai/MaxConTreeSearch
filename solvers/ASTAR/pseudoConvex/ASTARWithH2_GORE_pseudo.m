%A* tree search for pseudo-convex residuals
function [pk,vk,nodeNumber,uniqueNodeNumber, hInit, ubInit, levelMax,NOPruned, NODIBP] = ASTARWithH2_GORE_pseudo(A,b,c,d, x0, th, refine, BFSStart, runtimeLimit, newCheck)

if nargin < 7
    refine = 0;
end

if nargin < 8
    BFSStart = 0;
end

if nargin<9
    runtimeLimit = -1;
end


if nargin<10
    newCheck = 1;
end


dim1 = size(A,1)/numel(d);

acc = 0;
levelMax = 0;

NOPruned = 0;
NODIBP = 0;

uniqueNodeNumber = 1;
nodeNumber = 1;


[node.p, node.w, node.b] = minimaxPseudo(A,b,c,d,x0,th,dim1);
if (node.w<=acc)
    pk = node.p;
    vk = [];
    hInit = 0;
    ubInit = 0;
    return;
end


node.v = [];                    % Violation set.
N = numel(d);

[node.h, v, p] = computeLbAndUb_pseudo(A,b,c,d, [1:N], node.p, node.w, node.b, th, dim1, acc);     % Heuristic value.(X,y, H, v, x0, cyc, th)

%refine the initial upper bound
if refine == 1 && numel(v)-node.h > 5 && numel(x0) > 3
    [v,p] = IBCO_Refine_Quasi(A,b,c,d,p,th,node.h,v);
end

hInit = node.h;
U_cnt = numel(v);

pk = p;
vk = v;

disp(['ubInit = ' num2str(numel(vk)) ', lbInit = ' num2str(hInit)]);

ubInit = U_cnt;

if node.h >= U_cnt
    pk = p;
    vk = v;
    return;
end

node.f = node.h; % Evaluation value.

node.solved = 0;
minOUT = U_cnt;

q = node;

priof = node.f;
priow = node.w;
% List of checked bases.
dictS = java.util.Hashtable;

while 1
    
    [mn,~] = min(priof, [], 2);
    id = find(priof == mn);
    [~, mp] = min(priow(id));
    mp = id(mp);
    
    node = q(mp(1));
    if(node.solved)
        pk = node.p;
        vk = node.v;
        return;
    end
    q(mp) = [];
    
    priof(mp) = [];
    priow(mp) = [];
    
    if runtimeLimit > 0
        if toc > runtimeLimit
            disp(['reaching time limit of ' num2str(runtimeLimit)]);
            pk = p;
            vk = node.v;
            return;
        end
    end
    
    [childPoolTmp, bRest, bsPoolTmp, NOPruned, NODIBP,  q, priof, priow, dictS, pk, vk, minOUT, returnFlag] = GORE_Pseudo(A,b,c,d,dim1,acc,th, pk, vk, minOUT,node, q, priof, priow, dictS, NOPruned, NODIBP, N, newCheck);
    

    if returnFlag == 1
        return;
    end
    
    % Generate child nodes.
    for i=1:numel(bRest)
        
        child = childPoolTmp{bRest(i)};
        bs = bsPoolTmp{bRest(i)};
        nodeNumber = nodeNumber+1;
        
        H = 1:N;
        H(child.v) = [];
        
        uniqueNodeNumber = uniqueNodeNumber + 1;
        
        if BFSStart == 0
            [child.h, v, p] = computeLbAndUb_pseudo(A,b,c,d, H, child.p, child.w, bs, th, dim1, acc);     % Heuristic value.(X,y, H, v, x0, cyc, th)
        else
            hMax = floor((N-numel(child.v))/(numel(child.p)+1));
            if (numel(child.v)+hMax) >= minOUT                
                [child.h, v, p] = computeLbAndUb_pseudo(A,b,c,d, H, child.p, child.w, bs, th, dim1, acc);     % Heuristic value.(X,y, H, v, x0, cyc, th)

            else
                child.h = 0;
                p = child.p;
                
                
                HH = reshapeIdx(H,dim1);
                resTmp = resPseudoConvex(A(HH,:), b(:,H), c(:,H), d(H), p, th, dim1);
                v = find(abs(resTmp)>acc);
            end
        end

        UN_cnt = numel(v);
        child.h = max(child.h, node.h-1);
        child.f = numel(child.v)+child.h;
        
        if  minOUT > UN_cnt + numel(child.v)
            vk = [child.v, H(v)];
            pk = p;
            
            %refine the initial upper bound
            if refine == 1 && numel(vk)-child.f > 5 && numel(x0) > 3
               [vk,pk] = IBCO_Refine_Quasi(A,b,c,d,pk,th,child.f,vk);
            end


            
            disp(['after heuristic in genChild, vk -> ' num2str(numel(vk))]);
            minOUT = numel(vk);
            idx = priof > minOUT;
            
            q(idx) = [];
            priof(idx) = [];
            priow(idx) = [];
        end
        
        if numel(child.v) > levelMax
            levelMax = numel(child.v);
        end
        
        if child.f > minOUT
            continue;
        elseif child.h == UN_cnt
            % Reached a to the solution;
            child.solved = 1;
            child.p = p;
            child.w = th;
            child.v = sort([child.v, H(v)]);
        else
            child.solved = 0;
        end
        
        q = [ q, child ];
        priof = [priof, child.f];
        priow = [priow, child.w];
    end
end

end
