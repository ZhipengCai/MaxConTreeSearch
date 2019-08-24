function [ lwb, vioOut, solOut] = computeLbAndUb(X,y, H, solInit, fvalInit, idxActiveInit, th)

%compute the initial upper bound
Anew = X(H,:);
bnew = y(H);

resInit =  abs(Anew*solInit-bnew);
vioOut = find(resInit>th);

sol = solInit;
solOut = solInit;
ubOut = numel(vioOut);


%initialize the output
lwb = 0;

idxActive = idxActiveInit;
fval = fvalInit;

AOut = [];
bOut = [];


while(fval>th)
    
    AOut = [AOut; Anew(idxActive,:)];
    bOut = [bOut; bnew(idxActive)];
    
    Anew(idxActive,:) = [];
    bnew(idxActive) = [];
    
    resOut = abs(AOut(1:end-numel(idxActive),:)*sol-bOut(1:end-numel(idxActive)));
    idxMoveBack = find(resOut<=fval);
    
    
    Anew = [Anew; AOut(idxMoveBack,:)];
    bnew = [bnew; bOut(idxMoveBack)];
    AOut(idxMoveBack,:) = [];
    bOut(idxMoveBack) = [];
    
    if(numel(bnew) < numel(sol)+1)
        break;
    end
    
    %solve the constrained minimax problem (return the maximum residual and indices of the active constraint)
    [sol, fval, idxActive] = myFitTchebycheff(Anew,bnew, sol);
end


resOut = abs(X(H,:)*sol-y(H));
vioTmp = find(resOut>th);

if numel(vioTmp) < ubOut
    solOut = sol;
    vioOut = vioTmp;
    ubOut = numel(vioTmp);
end

j = 1;


dp1 = numel(sol)+1; %d+1
while j <= numel(bOut)
    if numel(bnew) < dp1
        NO2Add = min(dp1-numel(bnew), numel(bOut)-j+1);
        
        Anew = [Anew; AOut(j:(j+NO2Add-1),:)];
        bnew = [bnew; bOut(j:(j+NO2Add-1))];
        
        if numel(bnew) < dp1
            break;
        end
        
        j = j + NO2Add;
    else
        Anew = [Anew; AOut(j,:)];
        bnew = [bnew; bOut(j)];           
        j = j + 1;
    end
    
    
    %solve the constrained minimax problem (return the maximum residual and indices of the active constraint)
    [sol, fval, idxActive] = myFitTchebycheff(Anew,bnew, sol);
    
    if(fval>th)
        lwb = lwb+1;
        if ubOut<=lwb
            return;
        end
        
        Anew(idxActive,:) = [];
        bnew(idxActive) = [];
        
        
    end
end

resOut = abs(X(H,:)*sol-y(H));
vioTmp = find(resOut>th);

if numel(vioTmp) < ubOut
    solOut = sol;
    vioOut = vioTmp;
    ubOut = numel(vioTmp);
end

end
