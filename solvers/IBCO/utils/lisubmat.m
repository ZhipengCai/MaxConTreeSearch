function [Xsub,idx]=lisubmat(X,rowsOrCols,tol)
if ~nnz(X) %X has no non-zeros and hence no independent columns
    Xsub=[]; idx=[];
    return
end
if nargin<2
    rowsOrCols = 'cols'; tol=1e-10; 
elseif nargin<3
    tol=1e-10; 
end

if(strcmp(rowsOrCols, 'rows'))
    X=X';
    [Q, R, E] = qr(X,0);
    if ~isvector(R)
        diagr = abs(diag(R));
    else
        diagr = R(1);
    end
    %Rank estimation
    r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
    idx=sort(E(1:r));
    Xsub=X(:,idx)';
elseif (strcmp(rowsOrCols, 'cols'))
        
    [Q, R, E] = qr(X,0);
    if ~isvector(R)
        diagr = abs(diag(R));
    else
        diagr = R(1);
    end
    %Rank estimation
    r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
    idx=sort(E(1:r));
    Xsub=X(:,idx);

end

end