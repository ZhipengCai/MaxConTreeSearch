function [ lwb, vioOut, solOut] = computeLbAndUb_pseudo(A,b,c,d, H, solInit, fvalInit, idxActiveInit, th, dim1, acc)

HH = reshapeIdx(H,dim1);
%compute the initial upper bound
Anew = A(HH,:);
bnew = b(:,H);
cnew = c(:,H);
dnew = d(H);

resInit =  resPseudoConvex(Anew, bnew, cnew, dnew, solInit, th, dim1);

vioOut = find(resInit>acc);


sol = solInit;
solOut = solInit;
ubOut = numel(vioOut);


%initialize the output
lwb = 0;

idxActive = idxActiveInit;
fval = fvalInit;
AOut = [];
bOut = [];
cOut = [];
dOut = [];

while(fval>acc)
    
    idxActive2 = reshapeIdx(idxActive,dim1);
    
    AOut = [AOut; Anew(idxActive2,:)];
    bOut = [bOut, bnew(:,idxActive)];
    cOut = [cOut, cnew(:,idxActive)];
    dOut = [dOut, dnew(idxActive)];
    
    %remove active points
    Anew(idxActive2,:) = [];
    bnew(:,idxActive) = [];
    cnew(:,idxActive) = [];
    dnew(idxActive) = [];
        
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
    
    if(numel(dnew) < numel(sol)+1)
        break;
    end
    %solve the constrained minimax problem (return the maximum residual and indices of the active constraint)
    [sol, fval, idxActive] = minimaxPseudo(Anew,bnew,cnew,dnew, sol,th,dim1);

end


resOut = resPseudoConvex(A(HH,:), b(:,H), c(:,H), d(H), sol, th, dim1);

vioTmp = find(resOut>acc);

if numel(vioTmp) < ubOut
    solOut = sol;
    vioOut = vioTmp;
    ubOut = numel(vioTmp);
end

j = 1;

dp1 = numel(sol)+1; %d+1
while j <= numel(dOut)
    if numel(dnew) < dp1
        NO2Add = min(dp1-numel(dnew), numel(dOut)-j+1);
        
        idxTmp1 = j:(j+NO2Add-1);
        idxTmp2 = reshapeIdx(idxTmp1, dim1);
        
        Anew = [Anew; AOut(idxTmp2,:)];
        bnew = [bnew, bOut(:, idxTmp1)];
        cnew = [cnew, cOut(:, idxTmp1)];
        dnew = [dnew, dOut(idxTmp1)];
        
        if numel(dnew) < dp1
            break;
        end
        
        j = j + NO2Add;
    else
        idxTmp1 = j;
        idxTmp2 = reshapeIdx(idxTmp1, dim1);
        
        Anew = [Anew; AOut(idxTmp2,:)];
        bnew = [bnew, bOut(:, idxTmp1)];
        cnew = [cnew, cOut(:, idxTmp1)];
        dnew = [dnew, dOut(idxTmp1)];
        
        j = j + 1;
    end
        
        
    [sol, fval, idxActive] = minimaxPseudo(Anew,bnew,cnew,dnew, sol,th,dim1);
      
    
    
    if(fval>acc)
        lwb = lwb+1;
        if ubOut<=lwb
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
    vioOut = vioTmp;
end


end