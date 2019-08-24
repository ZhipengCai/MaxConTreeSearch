function [v,theta] = IBCO_Refine_Linear(X,y,theta0,th,lb,v)
ub = numel(v);
config.QThresh = 1e-9;
[A, b, c, d] = genQuasiconvexMatrixLinear(X, y);

theta = theta0;
N = numel(d);
d1 = numel(b)/N;
d2 = size(A,2);

%initialize sedumi model for constrains in socp problem for sx (only need to reset model.obj in BCO function)
epsilon_ac = th-config.QThresh; %tighten epsilon by accuracy to remove the effect of numerical issue
modelsx = genModelSocp_BCO_sedumi(A,b,c,d,epsilon_ac,N,d1,d2);


deltaMax = N-lb;
deltaMin = N-ub;
delta = floor((deltaMax+deltaMin)/2);

while deltaMax>(deltaMin+1)
    %step 1:
    [delta,deltaMin,deltaMax,theta] = BCO_v2(A, b, c, d, modelsx, theta, deltaMin, deltaMax, delta, config, th,N, d2, 'general');
end

if deltaMin > N-ub
    
    [~,~,inliers] = compute_residuals_l2(A,b,c,d,theta,th);
    if N-numel(inliers) < ub
        v = 1:N;
        v(inliers) = [];
    end
end
end