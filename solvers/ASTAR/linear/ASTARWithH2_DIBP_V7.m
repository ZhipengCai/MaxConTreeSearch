function [pk,vk,nodeNumber,uniqueNodeNumber, hInit, ubInit, levelMax, NOPruned, NODIBP] = ASTARWithH2_DIBP_V7(X,y, x0, th, refine, BFSStart, runtimeLimit)
if nargin < 5
    refine = 0;
end

if nargin < 6
    BFSStart = 0;
end

if nargin<7
    runtimeLimit = -1;
end

levelMax = 0;
NOPruned = 0;
NODIBP = 0;

N = size(X,1);

uniqueNodeNumber = 1;
nodeNumber = 1;

disp('initial Tcheb fit...');
[node.p, node.w, node.b] = myFitTchebycheff(X, y, x0);
if (node.w<=th)
    pk = node.p;
    vk = [];
    hInit = 0;
    ubInit = 0;
    return;
end


disp('initial heuristic...');
node.v = [];                    % Violation set.
[node.h, v, p] = computeLbAndUb(X,y, [1:N], node.p, node.w, node.b, th);   

%refine the initial upper bound
if refine == 1 && numel(v)-node.h > 5 && numel(x0) > 3
    [v,p] = IBCO_Refine_Linear(X,y,p,th,node.h,v);
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
    disp(['before return 1, vk -> ' num2str(numel(vk))]);
    return;
end

node.f = node.h; % Evaluation value.

node.solved = 0;
minOUT = U_cnt;

q = node;

priof = node.f;
priow = node.w;
dictS = java.util.Hashtable;
dictB = java.util.Hashtable;

while numel(q)>0
    if runtimeLimit > 0 
        if toc > runtimeLimit
            disp(['reaching time limit of ' num2str(runtimeLimit)]);
            return;
        end
    end
    
    [mn,~] = min(priof, [], 2);
    id = find(priof == mn);
    [~, mp] = min(priow(id));
    mp = id(mp);
    
    node = q(mp(1));
    if(node.solved)
        pk = node.p;
        vk = node.v;
        disp(['before return 2, vk -> ' num2str(numel(vk))]);
        return;
    end
    q(mp) = [];
    
    priof(mp) = [];
    priow(mp) = [];
    
    [childPoolTmp, bRest, bsPoolTmp, NOPruned, NODIBP,  q, priof, priow, dictS, dictB, pk, vk, minOUT, returnFlag] = DIBP_v7(X,y,th, pk, vk, minOUT,node, q, priof, priow, dictS, dictB, NOPruned, NODIBP, N);
    
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
        
            [child.h, v, p] = computeLbAndUb(X,y, H, child.p, child.w, bs, th);
        
        else
            hMax = floor((N-numel(child.v))/(numel(child.p)+1));
            if (numel(child.v)+hMax) >= minOUT
                [child.h, v, p] = computeLbAndUb(X,y, H, child.p, child.w, bs, th);
            else
                child.h = 0;
                p = child.p;
                resTmp = X(H,:)*p-y(H);
                v = find(abs(resTmp)>th);
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
                [vk,pk] = IBCO_Refine_Linear(X,y,pk,th,child.f,vk);
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
            disp(['reach level ' num2str(levelMax)]);
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
