%% GORE for A* branch pruning, lb and ub are refering to the number of outliers
function [lb, ubOut, solOut] = heuristicDIBP(X,y, H, idxEnforce, solInit, th, ub)
minSize = numel(solInit)+1-numel(idxEnforce);

%Initial solution
Anew = X(H,:);
bnew = y(H);

AEnforce = X(idxEnforce,:);
bEnforce = y(idxEnforce);

%solve the constrained minimax problem (return the maximum residual and indices of the active constraint)
[sol, fval, idxActive] = myFitTchebycheff_constrained(Anew, bnew, AEnforce, bEnforce, th, solInit);

%initialize the output
lb = 0;
solOut = solInit;
ubOut = ub;

AOut = [];
bOut = [];
while(fval>th)
    

    AOut = [AOut; Anew(idxActive,:)];
    bOut = [bOut; bnew(idxActive)];
    
    %remove active points
    Anew(idxActive,:) = [];
    bnew(idxActive) = [];
    
    resOut = abs(AOut(1:end-numel(idxActive),:)*sol-bOut(1:end-numel(idxActive)));
    idxMoveBack = find(resOut<=fval);
    Anew = [Anew; AOut(idxMoveBack,:)];
    bnew = [bnew; bOut(idxMoveBack)];
    AOut(idxMoveBack,:) = [];
    bOut(idxMoveBack) = [];
    
    if numel(bnew) < minSize
        break;
    end

[sol, fval, idxActive] = myFitTchebycheff_constrained(Anew, bnew, AEnforce, bEnforce, th, sol);

end


resOut = abs(X(H,:)*sol-y(H));
vioTmp = find(resOut>th);

if numel(vioTmp) < ubOut
    solOut = sol;
    ubOut = numel(vioTmp);
end

j = 1;

while j <= numel(bOut)
    if numel(bnew) < minSize
        NO2Add = min(minSize-numel(bnew), numel(bOut)-j+1);
        
        Anew = [Anew; AOut(j:(j+NO2Add-1),:)];
        bnew = [bnew; bOut(j:(j+NO2Add-1))];
        
        if numel(bnew) < minSize
            break;
        end
        
        j = j + NO2Add;
    else
        Anew = [Anew; AOut(j,:)];
        bnew = [bnew; bOut(j)];           
        j = j + 1;
    end
    
    [sol, fval, idxActive] = myFitTchebycheff_constrained(Anew, bnew, AEnforce, bEnforce, th, sol);

    
    if(fval>th)
        lb = lb+1;
        if ubOut<lb
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
    ubOut = numel(vioTmp);
end

end
