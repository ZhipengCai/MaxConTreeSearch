%***********************************************************************
%                     Demo program for paper:
%        "Consensus Maximization Tree Search Revisited"
%         Authors: Zhipeng Cai, Tat-Jun Chin, Vladlen Koltun 
%***********************************************************************

%-----------------------------------------------------------------------
%                           WARNING:
% - This demo is free for non-commercial academic use. Any commercial use is strictly prohibited without the authors' consent. 
% Please acknowledge the authors by citing:
% @article{cai2019consensus,
%  title={Consensus Maximization Tree Search Revisited},
%  author={Cai, Zhipeng and Chin, Tat-Jun and Koltun, Vladlen},
%  journal={arXiv preprint arXiv:1908.02021},
%  year={2019}
% }
%in any academic publications that have made use of this package or part of it.
%
%
%
% - This code was tested on a 64 bit Ubuntu 14.04 with MATLAB R2018b
% - To use this demo, just run function demo() in 'demo.m'
% - The time limit for all methods can be adjusted by changing "timeLimitFund"(for fundamental matrix estimation) and "timeLimitHomo" (for homography estimation).
%----------------------------------------------------------------------- 
clear;
close all;

%% important settings
%reduce the time limit to 20s/200s to make the demo finish quickly, u can adjust
%the time limit if u want.
timeLimitFund = 20; 
timeLimitHomo = 200;
selectedMethod = [5,4,3,2,1];

warning('off', 'all');
%% method list
methodList = {'A*', 'A*-TOD', 'A*-NAPA', 'A*-NAPA-TOD', 'A*-NAPA-DIBP'};

%% fundamental matrix estimation (the code for other experiments will be released after publication)
dataFolder = './Fund/data/';
dataList = {'KITTI_104_108', 'KITTI_198_201', 'KITTI_417_420', 'KITTI_579_582', 'KITTI_738_742'};

epsilonFund = 0.03; %inlier threshold

for i = 1:numel(dataList)
    data = char(dataList(i));
    dataFile = [dataFolder data '.mat'];
    disp(['reading: ' dataFile '...']);
    load(dataFile);
    %generate linearized data from sift matches
    [x, y] = genMatrixLinearizeFundamental(data.x1, data.x2);
    
    d = size(x,2);
    solInit = rand(d,1);
    %execute each method
    for j = 1:numel(selectedMethod)
        idxMethod = selectedMethod(j);
        disp('--------------------------------------------------------');
        disp(['executing ' methodList{idxMethod} '... N = ' num2str(numel(y)) '; epsilon = ' num2str(epsilonFund) ', runtimeLimit = ' num2str(timeLimitFund) 's']);
        [sol, outl, UniqueNodeNumber, hInit, ubInit, levelMax, NODIBP, runtime] = linearFitV2(x,y,solInit,epsilonFund,methodList{idxMethod},timeLimitFund);
        if idxMethod == 5
            outl_DIBP = outl;
        end
        disp([methodList{idxMethod} ' finished--runtime = ' num2str(runtime) ', NUN = ' num2str(UniqueNodeNumber), ', NOBP = ' num2str(NODIBP), ', levelMax = ' num2str(levelMax)]);
        disp('--------------------------------------------------------');
    end
    
    disp(['number of outliers o = ' num2str(numel(outl_DIBP))] );
    inls = 1:numel(y);
    inls(outl_DIBP) = [];
    plot_match(data.matches, [data.matches.X1; data.matches.X2], inls, 1, 1000);
    if i<numel(dataList)
        disp(['press any button to go to the next data']);
        pause;
    end
end


%% Homography estimation
dataFolder = './Homo/data/';
dataList = {'adam', 'city', 'Boston', 'Brussels', 'BruggeTower'};

epsilonHomoBasic = 4; %inlier threshold

for i = 1:numel(dataList)
    data = char(dataList(i));
    dataFile = [dataFolder data '.mat'];
    disp(['reading: ' dataFile '...']);
    load(dataFile);
    %generate linearized data from sift matches
    [A,b,c,d] = genMatrixHomography(data.x1, data.x2);

    d1 = size(A,2);
    solInit = rand(d1,1);
    
    epsilonHomo = epsilonHomoBasic*data.T2(1,1);
    
    %execute each method
    for j = 1:numel(selectedMethod)
        idxMethod = selectedMethod(j);
        disp('--------------------------------------------------------');
        disp(['executing ' methodList{idxMethod} '... N = ' num2str(numel(d)) '; epsilon = ' num2str(epsilonHomoBasic) ' pixels, runtimeLimit = ' num2str(timeLimitHomo) 's']);
        [sol, outl, UniqueNodeNumber, hInit, ubInit, levelMax, NODIBP, runtime] = pseudoConvFit(A,b,c,d,solInit,epsilonHomo,methodList{idxMethod},timeLimitHomo);
        if idxMethod == 5
            outl_DIBP = outl;
        end
        disp([methodList{idxMethod} ' finished--runtime = ' num2str(runtime) ', NUN = ' num2str(UniqueNodeNumber), ', NOBP = ' num2str(NODIBP), ', levelMax = ' num2str(levelMax)]);
        disp('--------------------------------------------------------');
    end
    
    disp(['number of outliers o = ' num2str(numel(outl_DIBP))] );
    inls = 1:numel(d);
    inls(outl_DIBP) = [];
    plot_match(data.matches, [data.matches.X1; data.matches.X2], inls, 1, 1000);
    if i<numel(dataList)
        disp(['press any button to go to the next data']);
        pause;
    end
end

