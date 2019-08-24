function [lb, ubOut, solOut] = heuristicDIBP_Pseudo(A,b,c,d, H, idxEnforce, solInit, th, dim1, acc, ub)

minSize = numel(solInit)+1-numel(idxEnforce);
HH = reshapeIdx(H,dim1);
Anew = A(HH,:);
bnew = b(:,H);
cnew = c(:,H);
dnew = d(H);


idxEnforce2 = reshapeIdx(idxEnforce,dim1);
AEnforce = A(idxEnforce2,:);
bEnforce = b(:,idxEnforce);
cEnforce = c(:,idxEnforce);
dEnforce = d(idxEnforce);


[res, sol, idxActive] = constrained_Minimax_Pseudo(Anew, bnew, cnew, dnew, AEnforce, bEnforce, cEnforce, dEnforce, th, dim1 , solInit, acc);


%initialize the output
lb = 0;
solOut = solInit;
ubOut = ub;

fval = res(idxActive(1));

AOut = [];
bOut = [];
cOut = [];
dOut = [];

AOutFix = [];
bOutFix = [];
cOutFix = [];
dOutFix = [];


while(fval>acc)
    
    idxActive2 = reshapeIdx(idxActive,dim1);
    
    if numel(idxActive) < minSize
        lb = lb+1;
        if ubOut<lb
            return;
        end
        moveBackFlag = 0;
        AOutFix = [AOutFix; Anew(idxActive2,:)];
        bOutFix = [bOutFix, bnew(:,idxActive)];
        cOutFix = [cOutFix, cnew(:,idxActive)];
        dOutFix = [dOutFix, dnew(idxActive)];
        
    else
        moveBackFlag = 1;
        AOut = [AOut; Anew(idxActive2,:)];
        bOut = [bOut, bnew(:,idxActive)];
        cOut = [cOut, cnew(:,idxActive)];
        dOut = [dOut, dnew(idxActive)];
        
    end
    
    Anew(idxActive2,:) = [];
    bnew(:,idxActive) = [];
    cnew(:,idxActive) = [];
    dnew(idxActive) = [];
    
    if moveBackFlag == 1
        resOut = resPseudoConvex(AOut(1:end-numel(idxActive2),:), bOut(:,1:end-numel(idxActive)), cOut(:, 1:end-numel(idxActive)), dOut(1:end-numel(idxActive)), sol, th, dim1);
        idxMoveBack = find(resOut<=fval);
        
        if numel(idxMoveBack)>0
            idxMoveBack2 = reshapeIdx(idxMoveBack, dim1);
            
            Anew = [Anew; AOut(idxMoveBack2,:)];
            bnew = [bnew, bOut(:,idxMoveBack)];
            cnew = [cnew, cOut(:,idxMoveBack)];
            dnew = [dnew, dOut(idxMoveBack)];
            
            AOut(idxMoveBack2,:) = [];
            bOut(:,idxMoveBack) = [];
            cOut(:,idxMoveBack) = [];
            dOut(idxMoveBack) = [];
        end
    end
    if numel(dnew) < minSize
        break;
    end
    
    
    [res, sol, idxActive] = constrained_Minimax_Pseudo(Anew, bnew, cnew, dnew, AEnforce, bEnforce, cEnforce, dEnforce, th,dim1, sol, acc);   
    
    fval = res(idxActive(1));
    
end


resOut = resPseudoConvex(A(HH,:), b(:,H), c(:,H), d(H), sol, th, dim1);
vioTmp = find(resOut>acc);

if numel(vioTmp) < ubOut
    solOut = sol;
    ubOut = numel(vioTmp);
end

j = 1;

while j <= numel(dOut)
    if numel(dnew) < minSize
        NO2Add = min(minSize-numel(dnew), numel(dOut)-j+1);
        
        idxTmp1 = j:(j+NO2Add-1);
        idxTmp2 = reshapeIdx(idxTmp1, dim1);
        
        Anew = [Anew; AOut(idxTmp2,:)];
        bnew = [bnew, bOut(:, idxTmp1)];
        cnew = [cnew, cOut(:,idxTmp1)];
        dnew = [dnew, dOut(idxTmp1)];
        
        if numel(dnew) < minSize
            break;
        end
        
        j = j + NO2Add;
    else
        
        idxTmp1 = j;
        idxTmp2 = reshapeIdx(idxTmp1, dim1);
        
        Anew = [Anew; AOut(idxTmp2,:)];
        bnew = [bnew, bOut(:, idxTmp1)];
        cnew = [cnew, cOut(:,idxTmp1)];
        dnew = [dnew, dOut(idxTmp1)];
        
        j = j + 1;
    end
    
    
    [res, sol, idxActive] = constrained_Minimax_Pseudo(Anew, bnew, cnew, dnew, AEnforce, bEnforce, cEnforce, dEnforce, th,dim1, sol,acc);
    
    
    fval = res(idxActive(1));
    
    
    if(fval>acc)
        lb = lb+1;
        if ubOut<lb
            return;
        end
        idxActive2 = reshapeIdx(idxActive,dim1);
        
        %remove active points
        Anew(idxActive2,:) = [];
        bnew(:,idxActive) = [];
        cnew(:,idxActive) = [];
        dnew(idxActive) = [];
        
    end
end

resOut = resPseudoConvex(A(HH,:), b(:,H), c(:,H), d(H), sol, th, dim1);
vioTmp = find(resOut>acc);

if numel(vioTmp) < ubOut
    solOut = sol;
    ubOut = numel(vioTmp);
end

end
