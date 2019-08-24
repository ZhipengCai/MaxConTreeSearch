function [theta, runtime] = IBCO_v2(A,b,c,d,theta0,th,config,problem)
if nargin < 8
    problem = 'general';
end

% From theta0, compute u0,s0, v0
%finds the number of inliers
theta = theta0;
N = numel(d);
d1 = numel(b)/N;
d2 = size(A,2);

%initialize sedumi model for constrains in socp problem for sx (only need to reset model.obj in BCO function)
epsilon_ac = th-config.QThresh; %tighten epsilon by accuracy to remove the effect of numerical issue
modelsx = genModelSocp_BCO_sedumi(A,b,c,d,epsilon_ac,N,d1,d2);

tic

%step 0:
nInls = inlierCountQuasi_l2(A,b,c,d,theta,th);
          
deltaMax = N;
deltaMin = nInls;

[delta,deltaMin,deltaMax, theta] = BCO_v2(A, b, c, d, modelsx, theta, deltaMin, deltaMax, N, config, th, N,d2, problem);


%Find Inliers and exit if problem solved for 0 outlier situation
if deltaMin == N
    runtime = toc;
    return;
end

while deltaMax>(deltaMin+1)
    %step 1:
    %execute BCO algorithm
    %disp(['delta =', num2str(delta), '; deltaMax = ', num2str(deltaMax), '; deltaMin = ', num2str(deltaMin)]);
    [delta,deltaMin,deltaMax,theta] = BCO_v2(A, b, c, d, modelsx, theta, deltaMin, deltaMax, delta, config, th,N, d2, problem);
end
runtime = toc;
end