%Dimension insensitive branch pruning for Pseudo-convex residuals
%% useAdaptive: if set to 0, S_DIBP will always be 1
function [childPoolTmp, bRest, bsPoolTmp, NOPruned, NODIBP, q, priof, priow, dictS, dictB, pk, vk, minOUT, returnFlag] = DIBP_Pseudo_V2(A,b,c,d,dim1,acc,th, pk, vk, minOUT,node, q, priof, priow, dictS, dictB, NOPruned, NODIBP, N, useAdaptive)
if nargin < 20
    useAdaptive = 0;
end

idxNodeB = reshapeIdx(node.b,dim1);
r = resPseudoConvex(A(idxNodeB,:), b(:,node.b), c(:,node.b), d(node.b), pk, th, dim1);
[~, idc] = sort(abs(r), 'descend');

childPoolTmp = []; 
bRest = []; 
bsPoolTmp = [];
bRemove = [];

returnFlag = 0;
dd = numel(pk);

ubTmp = minOUT - numel(node.v);
NTmp = N-numel(node.v);

if useAdaptive == 1
    S_DIBP = max(1,ceil(dd+1+(dd+1-NTmp)/(ubTmp-1))); 
else
    S_DIBP = 1;
end

NONonAdjacent = 0;

for i=1:length(node.b)
    %%%%%%%%%%%%%%%% check for repetition %%%%%%%%%%%%%
    child.v = sort([node.v node.b(idc(i))]);
    entry = dictS.get(sprintf('%4ld', child.v));
    if numel(entry)
        bRemove = [bRemove, idc(i)];
        continue;
    end
    
    
    H = 1:N;
    H(child.v) = [];
    HH = reshapeIdx(H,dim1);
    
    [child.p, child.w,bs] = minimaxPseudo(A(HH,:),b(:,H),c(:,H),d(H),node.p,th,dim1);
    if child.w<=acc
        pk = child.p;
        vk = child.v;
        returnFlag = 1;
        return;
    end
    
    child.b = sort(H(bs));
    entry = dictB.get(sprintf('%4ld', child.b));
    if numel(entry)
        bRemove = [bRemove, idc(i)];
        continue;
    end
    
    
    %% non-adjacent basis
    resRepCheck = resPseudoConvex(A,b,c,d,child.p,th,dim1);
    vioReal = find(resRepCheck>child.w);
    vioSize = numel(vioReal);
    
    if vioSize < numel(child.v)
        dictS.put(sprintf('%4ld', child.v), 1);
        % bug fixed, only use NAPA for repetition check
%         disp(size(child.v))
%         disp(size(vioReal'))
        child.v = sort(vioReal);
        entry = dictS.get(sprintf('%4ld', child.v));
        if numel(entry)
            NONonAdjacent = NONonAdjacent+1;
            continue;
        end
        dictS.put(sprintf('%4ld', child.v), 1);
%         NONonAdjacent = NONonAdjacent+1;
%         continue;
    end
    
    childPoolTmp{idc(i)} = child;
    bRest = [bRest, idc(i)];  %save the unrepeated points and use them as the enforced inliers
    bsPoolTmp{idc(i)} = bs;
end

%there is at least one element in bRest
if numel(bRest) == 0
    return;
end

if numel(bRemove) >= S_DIBP
    H = 1:N;
    H([node.v node.b(bRemove)]) = []; %remove enforced inliers and the previously removed violation set
    NODIBP = NODIBP + 1;
    
    [h, ubDIBP, solDIBP] = heuristicDIBP_Pseudo(A, b, c, d, H, node.b(bRemove), pk, th, dim1, acc, ubTmp);
    
    if ubDIBP < ubTmp
        pk = solDIBP;
        resTmp = resPseudoConvex(A,b,c,d,solDIBP,th,dim1);
        vkTmp = find(resTmp > acc);
        if minOUT > vkTmp
            minOUT = numel(vkTmp);
            vk = vkTmp;
            disp(['during DIBP, vk -> ' num2str(numel(vk))]);
            
            ubTmp = minOUT - numel(node.v);
            
            if minOUT <= node.f
                returnFlag = 1;
                return;
            end
            
            idx = priof > minOUT;
            
            q(idx) = [];
            priof(idx) = [];
            priow(idx) = [];
        end
    end
    
    if h > minOUT - numel(node.v)
        bRest = [];
        NOPruned = NOPruned + numel(bRest)-i;
        return;
    end
    
    %% perform DIBP
    for i = 1:numel(bRest)-1
        
        dictS.put(sprintf('%4ld', childPoolTmp{bRest(i)}.v), 1);
        dictB.put(sprintf('%4ld', childPoolTmp{bRest(i)}.b), 1);
        
        H = 1:N;
        H([node.v node.b(bRest(1:i)) node.b(bRemove)]) = []; 
        NODIBP = NODIBP + 1;
        
        [h, ubDIBP, solDIBP] = heuristicDIBP_Pseudo(A, b, c, d, H, node.b([bRest(1:i), bRemove]), pk, th, dim1, acc, ubTmp);
        
        if ubDIBP < ubTmp
            pk = solDIBP;
            resTmp = resPseudoConvex(A,b,c,d,solDIBP,th,dim1);
            vkTmp = find(resTmp > acc);
            if minOUT > vkTmp
                minOUT = numel(vkTmp);
                vk = vkTmp;
                disp(['during DIBP, vk -> ' num2str(numel(vk))]);
                
                ubTmp = minOUT - numel(node.v);
                
                if minOUT <= node.f
                    returnFlag = 1;
                    return;
                end
                
                idx = priof > minOUT;
                
                q(idx) = [];
                priof(idx) = [];
                priow(idx) = [];
            end
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
    for i = 1:S_DIBP-numel(bRemove)-1
        dictS.put(sprintf('%4ld', childPoolTmp{bRest(i)}.v), 1);
        dictB.put(sprintf('%4ld', childPoolTmp{bRest(i)}.b), 1);
    end
    
    %% perform DIBP
    for i = S_DIBP-numel(bRemove):numel(bRest)-1
        
        
        dictS.put(sprintf('%4ld', childPoolTmp{bRest(i)}.v), 1);
        dictB.put(sprintf('%4ld', childPoolTmp{bRest(i)}.b), 1);
        
        H = 1:N;
        H([node.v node.b(bRest(1:i)) node.b(bRemove)]) = []; 
        NODIBP = NODIBP + 1;
        
        [h, ubDIBP, solDIBP] = heuristicDIBP_Pseudo(A, b, c, d, H, node.b([bRest(1:i), bRemove]), pk, th, dim1, acc, ubTmp);
        
        if ubDIBP < ubTmp
            pk = solDIBP;
            resTmp = resPseudoConvex(A,b,c,d,solDIBP,th,dim1);
            
            vkTmp = find(resTmp > acc);
            if minOUT > vkTmp
                minOUT = numel(vkTmp);
                vk = vkTmp;
                disp(['during DIBP, vk -> ' num2str(numel(vk))]);
                
                ubTmp = minOUT - numel(node.v);
                
                if minOUT <= node.f
                    returnFlag = 1;
                    return;
                end
                
                idx = priof > minOUT;
                
                q(idx) = [];
                priof(idx) = [];
                priow(idx) = [];
            end
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