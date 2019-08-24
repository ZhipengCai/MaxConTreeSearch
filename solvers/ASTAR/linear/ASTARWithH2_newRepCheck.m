function [pk,vk,nodeNumber,uniqueNodeNumber, hInit, ubInit, levelMax] = ASTARWithH2_newRepCheck(X,y, x0, th, refine, BFSStart, runtimeLimit, newCheck)

if nargin < 5
    refine = 0;
end

if nargin < 6
    BFSStart = 0;
end

if nargin<7
    runtimeLimit = -1;
end

if nargin<8
    newCheck = 1;
end

levelMax = 0;

N = size(X,1);

uniqueNodeNumber = 1;
nodeNumber = 1;

[node.p, node.w, node.b] = myFitTchebycheff(X, y, x0);
if (node.w<=th)
    pk = node.p;
    vk = [];
    hInit = 0;
    ubInit = 0;
    return;
end

node.v = [];                    % Violation set.
[node.h, v, p] = computeLbAndUb(X,y, [1:N], node.p, node.w, node.b, th);     % Heuristic value.(X,y, H, v, x0, cyc, th)

%refine the initial upper bound
if refine == 1 && numel(v)-node.h > 5 && numel(x0) > 3
    [v,p] = IBCO_Refine_Linear(X,y,p,th,node.h,v);
end

hInit = node.h;
U_cnt = numel(v);

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
    
    % Generate child nodes.
    for i=1:length(node.b)
        child.v = sort([node.v node.b(i)]);
        entry = dictS.get(sprintf('%4ld', child.v));
        if numel(entry)
            continue;
        end
        
        dictS.put(sprintf('%4ld', child.v), 1);
        
        nodeNumber = nodeNumber+1;
        
        H = 1:N;
        H(child.v) = [];
        [child.p, child.w,bs] = myFitTchebycheff(X(H,:),y(H), node.p);
        
        
        if child.w<=th
            pk = child.p;
            vk = child.v;
            return;
        end
        
        child.b = sort(H(bs));
        
        if newCheck == 1
            %% new repetition check
            vioReal = find(abs(X*child.p-y)>child.w);
            vioSize = numel(vioReal);
            
            if vioSize < numel(child.v)
                continue;
            end
        end
        
        
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
            minOUT = UN_cnt + numel(child.v);
            
            vk = [child.v, H(v)];
            pk = p;
            
            if refine == 1 && numel(vk)-child.f > 5 && numel(x0) > 3
                [vk,pk] = IBCO_Refine_Linear(X,y,pk,th,child.f,vk);
                minOUT = numel(vk);
            end
            
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

