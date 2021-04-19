function [childPoolTmp, bRest, bsPoolTmp, NOPruned, NODIBP, q, priof, priow, dictS, pk, vk, minOUT, returnFlag] = GORE_v5(X,y,th, pk, vk, minOUT,node, q, priof, priow, dictS, NOPruned, NODIBP, N, newCheck)
r = X(node.b, :)*pk - y(node.b);
[~, idc] = sort(abs(r), 'descend');

childPoolTmp = [];
bRest = []; 
bsPoolTmp = [];

returnFlag = 0;
d = numel(pk);

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
    
    [child.p, child.w,bs] = myFitTchebycheff(X(H,:), y(H), node.p);
    
    if child.w<=th
        % Reached a goal state.
        pk = child.p;
        vk = child.v;
        returnFlag = 1;
        return;
    end
    
    child.b = sort(H(bs));
    
    if newCheck == 1
        %% new repetition check
        vioReal = find(abs(X*child.p-y)>child.w);
        vioSize = numel(vioReal);
        
        if vioSize < numel(child.v)
            dictS.put(sprintf('%4ld', child.v), 1);
%             continue;
            % bug removed, only use NAPA for repetition check
            child.v = sort(vioReal');
            entry = dictS.get(sprintf('%4ld', child.v));
            if numel(entry)
                continue;
            end
            dictS.put(sprintf('%4ld', child.v), 1);
        end
    end
    
    %%%%%%%%%%%%%%%% repetition check ends, and results reused later %%%%%%%%%%%%%
    childPoolTmp{idc(i)} = child;
    bRest = [bRest, idc(i)];  %save the unrepeated points and use them as the enforced inliers
    bsPoolTmp{idc(i)} = bs;
end

if numel(bRest) < 2
    return;
else
    for i = 1:numel(bRest)-1
        dictS.put(sprintf('%4ld', childPoolTmp{bRest(i)}.v), 1);
        H = 1:N;
        H([node.v node.b(bRest(i))]) = []; %remove enforced inliers and the previously removed violation set
        NODIBP = NODIBP + 1;
        [h, ubDIBP, solDIBP] = heuristicDIBP(X, y, H, node.b(bRest(i)), pk, th, ubTmp);
        
        if ubDIBP < ubTmp
            pk = solDIBP;
            resTmp = abs(X*solDIBP-y);
            vk = find(resTmp > th+1e-8);
            minOUT = min(minOUT, numel(vk));
            
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
        
        if h > minOUT - numel(node.v)
            bRest = bRest(i);
            NOPruned = NOPruned + numel(bRest)-1;
            return;
        end
    end
    
    dictS.put(sprintf('%4ld', childPoolTmp{bRest(end)}.v), 1);
end

end