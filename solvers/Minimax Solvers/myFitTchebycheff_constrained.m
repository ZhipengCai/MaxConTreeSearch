function [xn,d,M] = myFitTchebycheff_constrained(X, Y, XCon, YCon, th, xn)

if nargin < 6
    xn = linsolve(XCon, YCon);
end


%number of parameters
n = size(X, 2);

%number of points
l = numel(Y);

%replace each constraint with an absolute value by two linear constraints
X = [X; -X];
Y = [Y; -Y];

rCon = abs(XCon*xn-YCon);
if sum(rCon>th)>=1
    xn = linsolve(XCon, YCon);
end

%same for the constraints
XCon = [XCon; -XCon];
YCon = [YCon; -YCon];


rCon = XCon*xn-YCon;


r = X*xn - Y;


%maximum value of the linear residuals
d = max(r, [], 1);
%find the indices (of the linear constraints) with current maximum
%residuals
M = find(abs(r-d)<1e-9);
%set the residual of those points to the maximum value
r(M) = d;

MCon = find(abs(rCon-th)<1e-9);
rCon(MCon) = th;

while(numel(M)+numel(MCon) < n+1)
    if d < 1e-9
        M = [];
        d = 0;
        return;
    end
    % compute the descent direction such that all maximal residuals
    % decrease in an equal speed
    A = [X(M, :); XCon(MCon,:)];
    c = -(A*A')\[ones(numel(M), 1); zeros(numel(MCon),1)];    
    y = A'*c;
        
    %compute the step size along the descent direction such that another
    %residual becomes maximal
    lambda = (d - r)./(X*y+1);
    lambda(isnan(lambda)) = 0;
    lambda(lambda<=0) = Inf;
    
    denom = XCon*y;
    lambdaCon = (th-rCon)./(denom);
    lambdaCon(abs(denom)<1e-8) = 0;
    lambdaCon(isnan(lambdaCon)) = 0;
    lambdaCon(lambdaCon<=0) = Inf;
    
    [lambda_j, j] = min(lambda);
    [lambdaCon_j, jCon] = min(lambdaCon);
    
    if lambda_j < lambdaCon_j
        xn = xn + lambda_j*y;
        
        M = [M; j];
        r = X*xn - Y;
        d = max(r, [], 1);
        r(M) = d;
        
        rCon = XCon*xn-YCon;
        rCon(MCon) = th;
    elseif lambda_j == lambdaCon_j
        xn = xn + lambda_j*y;
        
        M = [M; j];
        r = X*xn - Y;
        d = max(r, [], 1);
        r(M) = d;
        
        MCon = [MCon; jCon];
        rCon = XCon*xn-YCon;
        rCon(MCon) = th;
    else
        xn = xn + lambdaCon_j*y;
        
        r = X*xn - Y;
        d = max(r, [], 1);
        r(M) = d;
        
        MCon = [MCon; jCon];
        rCon = XCon*xn-YCon;
        rCon(MCon) = th;
    end
end

if n+1~=numel(M)+numel(MCon)
    [res, xn, M] = constrained_Minimax_Linear(X(1:size(X,1)/2,:), Y(1:size(X,1)/2,:), XCon(1:size(XCon,1)/2,:), YCon(1:size(XCon,1)/2,:), th,xn);
    d = res(M(1));
    return;
end

C = inv([[ones(numel(M), 1); zeros(numel(MCon),1)], [X(M', :); XCon(MCon',:)]]);
C1 = C(1, :);

itrn = 1e6;


while(sum(C1>=-1e-9) < n+1) && itrn
    if d < 1e-9
        M = [];
        d = 0;
        return;
    end
    
        if 1e6-itrn+1 > 100
            disp(['iter ' num2str(1e6-itrn+1) ' of stage 2...']);
            
        end
    
    [~, p1] = min(C1); %% Error Occuring ingeneral in p1-th position
    
    if abs(C(1, p1))<1e-9
        
        y = C(2:end, p1)/(-1e-9);
        
    else
        
        y = C(2:end, p1)/C(1, p1);
    
    end
    
    t = (d - r)./(X*y+1);
    t(isnan(t)) = 0;
    t(t<=0) = Inf;
    [t_j, j] = min(t);
    
    
    denom = XCon*y;
    tCon = (th-rCon)./denom;
    tCon(abs(denom)<1e-8) = 0;

    tCon(isnan(tCon)) = 0;
    tCon(tCon<=0) = Inf;
    [tCon_j, jCon] = min(tCon);
        
    if t_j < tCon_j
        xn = xn + t_j*y;
        
        lambda = [1, X(j, :)]*C;
        
        beta = p1;
    
        if beta<=numel(M)
            M(beta) = j;          %% In some cases beta replaced by p1 solves
        
            Cb = C(:, beta)/lambda(beta);
            C = C - (ones(n+1, 1)*lambda) .* (Cb*ones(1, n+1));
            C(:, beta) = Cb;
        else
            
            MCon(beta-numel(M)) = [];
            M(end+1) = j;
            
            Cb = C(:, beta)/lambda(beta);
            C = C - (ones(n+1, 1)*lambda) .* (Cb*ones(1, n+1));
            C(:, beta) = [];
            C = [C(:,1:numel(M)-1), Cb, C(:,numel(M):end)];
        end
        
        
        
        r = X*xn - Y;
        d = max(r);
        r(M) = d;
        
        rCon = XCon*xn-YCon;
        
        rCon(MCon) = th;
        C1 = C(1, :);
        
    else
        xn = xn + tCon_j*y;
        lambda = [0, XCon(jCon, :)]*C;
        
        beta = p1; % This huristic
        
        if beta>numel(M)
            MCon(beta-numel(M)) = [jCon];
            
            Cb = C(:, beta)/lambda(beta);
            C = C - (ones(n+1, 1)*lambda) .* (Cb*ones(1, n+1));
            C(:, beta) = Cb;
        else
            M(beta) = [];
            MCon = [MCon; jCon];
            
            Cb = C(:, beta)/lambda(beta);
            C = C - (ones(n+1, 1)*lambda) .* (Cb*ones(1, n+1));
            C = [C(:,1:beta-1), C(:,beta+1:end), Cb];
        end
        
        
      
        
        r = X*xn - Y;
        d = max(r);
        r(M) = d;
        
        rCon = XCon*xn-YCon;
        rCon(MCon) = th;
        
        C1 = C(1, :);
    end
    itrn = itrn - 1;
end

if d < 1e-9
    M = [];
    d = 0;
    return;
end

id = M > l;
M(id) = M(id) - l;
