%biconvex optimization
function [delta,deltaMin,deltaMax,theta0] = BCO_v2(A, b, c, d, modelsx, theta0, deltaMin, deltaMax, delta, config,epsilon,n,d2,problem)
if delta == n %solve socp when delta = 0
    y = ones(n,1);
    pars.fid=0;
    bt = [zeros(d2,1); -ones(n,1)];
    [~,sx1,~] = sedumi(modelsx.At,bt,modelsx.ct,modelsx.K,pars);
    theta = sx1(1:d2);
    currF = sum(sx1(d2+1:d2+n));
    %enforce the constraint on theta for fundamental matrix estimation
    if(strcmp(problem,'fun'))
        F = [theta(1:3)';theta(4:6)';[theta(7),1,theta(8)]];
        FPrime = enforceFundamentalConstraint(F);
        theta = [FPrime(1,:) FPrime(2,:) FPrime(3,1) FPrime(3,3)]';
        [~,sx1((d2+1):(d2+n)),~] = compute_residuals_l2(A,b,c,d,theta,epsilon);
        FProj = sum(sx1(d2+1:d2+n));
    else
        FProj = currF;
    end
    nInls = inlierCountQuasi_l2(A,b,c,d,theta,epsilon);
    if deltaMin<nInls
        theta0 = theta;
        deltaMin = nInls;
    end
    delta = deltaMin+floor((deltaMax-deltaMin)/3*2);
else %use biconvex optimization technique
    maxiter = 1e6;%safe guard to prevent infinite iterations
    currF = 1e9;
    theta = theta0;
    for i = 1:maxiter
        preF = currF;
        %1. fix sx, solve u
        [~,s,~] = compute_residuals_l2(A,b,c,d,theta,epsilon);
        y = solveY_BCO(s,delta);          % u can be solved in close form
        pars.fid=0;
        bt = [zeros(d2,1); -y];
        [~,sx1,info] = sedumi(modelsx.At,bt,modelsx.ct,modelsx.K,pars);
        
        theta = sx1(1:d2);
        [~,s,inliers] = compute_residuals_l2(A,b,c,d,theta,epsilon);
        currF = y'*s;
        
        if(i>1 && (abs(preF-currF)<1e-9 || preF<currF))
            if(strcmp(problem,'fun'))
                F = [theta(1:3)';theta(4:6)';[theta(7),1,theta(8)]];
                FPrime = enforceFundamentalConstraint(F);
                theta = [FPrime(1,:) FPrime(2,:) FPrime(3,1) FPrime(3,3)]';
                [~,s,inliers] = compute_residuals_l2(A,b,c,d,theta,epsilon);
                y = solveY_BCO(s,delta);          % u can be solved in close form
                FProj = y'*s;
            else
                FProj = currF;
            end
            break;
        end
    end
    
    nInls = numel(inliers);
    if deltaMin<nInls
        theta0 = theta;
        deltaMin = nInls;
    end
    if FProj > config.QThresh
         deltaMax = delta; 
    end
    delta = floor((deltaMax+deltaMin)/2);
end
end