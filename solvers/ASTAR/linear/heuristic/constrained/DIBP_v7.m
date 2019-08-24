%Dimension insensitive branch pruning with adaptive selection
%of the initial number of enforced inliers
function [childPoolTmp, bRest, bsPoolTmp, NOPruned, NODIBP, q, priof, priow, dictS,dictB, pk, vk, minOUT, returnFlag] = DIBP_v7(X,y,th, pk, vk, minOUT,node, q, priof, priow, dictS,dictB, NOPruned, NODIBP, N)
r = X(node.b, :)*pk - y(node.b);
[~, idc] = sort(abs(r), 'descend');

childPoolTmp = []; 
bRest = []; 
bsPoolTmp = [];
bRemoved = []; 
returnFlag = 0;
d = numel(pk);

ubTmp = minOUT - numel(node.v);
NTmp = N-numel(node.v);
S_DIBP = max(1,ceil(d+1+(d+1-NTmp)/(ubTmp-1))); 

NONonAdjacent = 0;

% find and remove all the repeated children
for i=1:length(node.b)
    %%%%%%%%%%%%%%%% check for repetition %%%%%%%%%%%%%
    child.v = sort([node.v node.b(idc(i))]);
    entry = dictS.get(sprintf('%4ld', child.v));
    if numel(entry)
        bRemoved = [bRemoved, idc(i)];
        continue;
    end
    
    
    H = 1:N;
    H(child.v) = [];
    
    [child.p, child.w,bs] = myFitTchebycheff(X(H,:), y(H), node.p);
    
    if child.w<=th
        % Reached a goal state.
        pk = child.p;
        vk = child.v;
        disp(['before return GORE, vk -> ' num2str(numel(vk))]);
        returnFlag = 1;
        return;
    end
    
    child.b = sort(H(bs));
    %repetition in different levels
    entry = dictB.get(sprintf('%4ld', child.b));
    if numel(entry)
        dictS.put(sprintf('%4ld', child.v), 1);
        bRemoved = [bRemoved, idc(i)];
        continue;
    end
    
    %% non-adjacent basis
    vioReal = find(abs(X*child.p-y)>child.w);
    vioSize = numel(vioReal);
    if vioSize < numel(child.v)
        dictS.put(sprintf('%4ld', child.v), 1);
        NONonAdjacent = NONonAdjacent+1;
        continue;
    end
    
    %%%%%%%%%%%%%%%% repetition check ends, and results reused later %%%%%%%%%%%%%
    childPoolTmp{idc(i)} = child;
    bRest = [bRest, idc(i)]; 
    bsPoolTmp{idc(i)} = bs;
end

if numel(bRemoved)+NONonAdjacent>d
    return;
end

%try to use bRemoved to prune
if numel(bRemoved) >= S_DIBP
    H = 1:N;
    H([node.v node.b(bRemoved)]) = []; %remove enforced inliers and the previously removed violation set
    NODIBP = NODIBP + 1;
    [h, ubDIBP, solDIBP] = heuristicDIBP(X, y, H, node.b(bRemoved), pk, th, ubTmp);
    if ubDIBP < ubTmp
        pk = solDIBP;
        resTmp = abs(X*solDIBP-y);
        vk = find(resTmp > th+1e-8);
        minOUT = min(minOUT, numel(vk));
        
        ubTmp = minOUT - numel(node.v);
        disp(['during GORE, vk -> ' num2str(numel(vk))]);
        
        if minOUT <= node.f
            returnFlag = 1;
            return;
        end
        
        idx = priof > minOUT;
        
        q(idx) = [];
        priof(idx) = [];
        priow(idx) = [];
    end
    
    if h > minOUT - numel(node.v)
        bRest = [];
        NOPruned = NOPruned + numel(bRest);
        return;
    end
    
    %     disp('case 1 pruning stage 2');
    for i = 1:numel(bRest)-1
        
        dictS.put(sprintf('%4ld', childPoolTmp{bRest(i)}.v), 1);
        dictB.put(sprintf('%4ld', childPoolTmp{bRest(i)}.b), 1);
        H = 1:N;
        H([node.v node.b(bRest(1:i)) node.b(bRemoved)]) = []; %remove enforced inliers and the previously removed violation set
        NODIBP = NODIBP + 1;
        [h, ubDIBP, solDIBP] = heuristicDIBP(X, y, H, node.b([bRest(1:i), bRemoved]), pk, th, ubTmp);
        if ubDIBP < ubTmp
            pk = solDIBP;
            resTmp = abs(X*solDIBP-y);
            vk = find(resTmp > th+1e-8);
            minOUT = min(minOUT, numel(vk));
            
            ubTmp = minOUT - numel(node.v);
            disp(['during GORE, vk -> ' num2str(numel(vk))]);
            
            if minOUT <= node.f
                returnFlag = 1;
                return;
            end
            
            idx = priof > minOUT;
            
            q(idx) = [];
            priof(idx) = [];
            priow(idx) = [];
        end
        
        if h > minOUT - numel(node.v)
            bRest = bRest(1:i);
            NOPruned = NOPruned + numel(bRest)-i;
            return;
        end
    end
    
    dictS.put(sprintf('%4ld', childPoolTmp{bRest(end)}.v), 1);
    dictB.put(sprintf('%4ld', childPoolTmp{bRest(end)}.b), 1);
    
else
    
    if S_DIBP-numel(bRemoved)-1>numel(bRest)
        for i = 1:numel(bRest)
            dictS.put(sprintf('%4ld', childPoolTmp{bRest(i)}.v), 1);
            dictB.put(sprintf('%4ld', childPoolTmp{bRest(i)}.b), 1);
        end
        return;
    end
    
    for i = 1:S_DIBP-numel(bRemoved)-1
        dictS.put(sprintf('%4ld', childPoolTmp{bRest(i)}.v), 1);
        dictB.put(sprintf('%4ld', childPoolTmp{bRest(i)}.b), 1);
    end
    
    %% perform DIBP
    for i = S_DIBP-numel(bRemoved):numel(bRest)-1
        dictS.put(sprintf('%4ld', childPoolTmp{bRest(i)}.v), 1);
        dictB.put(sprintf('%4ld', childPoolTmp{bRest(i)}.b), 1);
        
        H = 1:N;
        H([node.v node.b(bRest(1:i)) node.b(bRemoved)]) = []; %remove enforced inliers and the previously removed violation set
        NODIBP = NODIBP + 1;
        [h, ubDIBP, solDIBP] = heuristicDIBP(X, y, H, node.b([bRest(1:i), bRemoved]), pk, th, ubTmp);
        if ubDIBP < ubTmp
            pk = solDIBP;
            resTmp = abs(X*solDIBP-y);
            vk = find(resTmp > th+1e-8);
            minOUT = min(minOUT, numel(vk));
            
            ubTmp = minOUT - numel(node.v);
            disp(['during GORE, vk -> ' num2str(numel(vk))]);
            
            if minOUT <= node.f
                returnFlag = 1;
                return;
            end
            
            idx = priof > minOUT;
            
            q(idx) = [];
            priof(idx) = [];
            priow(idx) = [];
        end
        
        if h > minOUT - numel(node.v)
            bRest = bRest(1:i);
            NOPruned = NOPruned + numel(bRest)-i;
            return;
        end
    end
    
    dictS.put(sprintf('%4ld', childPoolTmp{bRest(end)}.v), 1);
    dictB.put(sprintf('%4ld', childPoolTmp{bRest(end)}.b), 1);
end



end