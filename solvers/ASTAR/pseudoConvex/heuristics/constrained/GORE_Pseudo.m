function [childPoolTmp, bRest, bsPoolTmp, NOPruned, NOGORE, q, priof, priow, dictS, pk, vk, minOUT, returnFlag] = GORE_Pseudo(A,b,c,d,dim1,acc,th, pk, vk, minOUT,node, q, priof, priow, dictS, NOPruned, NOGORE, N, newCheck)

idxNodeB = reshapeIdx(node.b,dim1);
r = resPseudoConvex(A(idxNodeB,:), b(:,node.b), c(:,node.b), d(node.b), pk, th, dim1);
[~, idc] = sort(abs(r), 'descend');

childPoolTmp = []; 
bRest = []; 
bsPoolTmp = [];

returnFlag = 0; 
dd = numel(pk);

ubTmp = minOUT - numel(node.v);
NTmp = N-numel(node.v);


% find and remove all the repeated children
for i=1:length(node.b)
    %%%%%%%%%%%%%%%% check for repetition %%%%%%%%%%%%%
    child.v = sort([node.v node.b(idc(i))]);
    entry = dictS.get(sprintf('%4ld', child.v));
    if numel(entry)
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
    if newCheck == 1
        %% new repetition check
        resRepCheck = resPseudoConvex(A,b,c,d,child.p,th,dim1);
        vioReal = find(resRepCheck>child.w);
        vioSize = numel(vioReal);
        
        if vioSize < numel(child.v)
            dictS.put(sprintf('%4ld', child.v), 1);
            % bug fixed, only use NAPA for repetition check
            child.v = sort(vioReal);
            entry = dictS.get(sprintf('%4ld', child.v));
            if numel(entry)
                continue;
            end
            dictS.put(sprintf('%4ld', child.v), 1);
%             continue;
        end
    end
    
    childPoolTmp{idc(i)} = child;
    bRest = [bRest, idc(i)]; 
    bsPoolTmp{idc(i)} = bs;
end

if numel(bRest) < 2
    return;
else
    for i = 1:numel(bRest)
        dictS.put(sprintf('%4ld', childPoolTmp{bRest(i)}.v), 1);
        H = 1:N;
        H([node.v node.b(bRest(i))]) = []; 
        NOGORE = NOGORE + 1;

        [h, ubDIBP, solDIBP] = heuristicDIBP_Pseudo(A, b, c, d, H, node.b(bRest(i)), pk, th, dim1, acc, ubTmp);

        if ubDIBP < ubTmp
            pk = solDIBP;
            resTmp = resPseudoConvex(A,b,c,d,solDIBP,th,dim1);
            vk = find(resTmp > acc);
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
            bRest = bRest(i);
            NOPruned = NOPruned + numel(bRest)-i;
            return;
        end
    end
    
    dictS.put(sprintf('%4ld', childPoolTmp{bRest(end)}.v), 1);
end

end