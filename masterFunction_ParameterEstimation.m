clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics and fit to experimental data.
%% Two approaches can be implemeneted: Genetic Algorithms and Hooke and Jeeves Algorithm.
%% The code reproduces figure 4 and results given in table of parameters.
%% All results are stored in Results_Fits.

folderName = horzcat('Results_Fits4'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end

%% Deffining the optimization method and parameters.
runningOptimization =    0; % 0 plots the results with the saved best fit. 1 runs the optimization process.
optimizationMethod = 'GA'; % Options are 'GA' for Genetic Algorithms or 'PS' for Hooke and Jeeves Algorithm.
PopulationSize = 20; generations = 20;

if runningOptimization ==1
    
    %% Optmization implementation
    % To facilitate the optimization process the parameters are scaled to
    % integer values. Elongation rates are scaled as ke_scale = ke_real *10
    % Initiation rates are scaled as ki_scaled = ki_real*1e4;
    ke_min_1 = 80;  ke_max_1 =110;
    ki_min_1 = 10;  ki_max_1 =80;
    ki_min_2 = 10;  ki_max_2 =80;
    ki_min_3 = 10;  ki_max_3 =90;
%     initialMatrix = [98.0000   23.0000   36.5000   53.5000];
    initialMatrix = [106    22    50    66];

    %% loading experimental intensities.
    [experimentalData] = experimentalIntesityDistribution;
    % gene =1; % KDM5B
    % gene =2; % b-act
    % gene =3; % H2B
    
    %% deffining parameter ranges
    lb = [ke_min_1,  ki_min_1, ki_min_2, ki_min_3];
    ub = [ke_max_1,  ki_max_1, ki_max_2, ki_max_3];
    nvars= length(lb);
    plottingCondition =0;
    IntCon = [1,2,3,4]; % parameters with integer values.
    %% Parameter matrix
    fitFnc = @(x)fitData(x,experimentalData, plottingCondition,folderName);
    %% Running the Pattern Search Algorithm
    switch optimizationMethod
        case 'PS'
            options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf,'MaxIterations',1000,'MaxTime',3600*2);%
            x = patternsearch(fitFnc,initialMatrix,[],[],[],[],lb,ub,options);
        case 'GA'
            optionsGA = optimoptions('ga','PlotFcn',{@gaplotbestf},'MaxGenerations',generations,'PopulationSize',PopulationSize,'ConstraintTolerance',1e-6,'InitialPopulationMatrix',initialMatrix,'Display', 'iter'); %
            [x,fval] = ga(fitFnc,nvars,[],[],[],[],lb,ub,[],IntCon,optionsGA);
    end
    
    %% Reporting Results
    plottingCondition=1;
    [experimentalData] = experimentalIntesityDistribution ;
    fitVal = fitData(x,experimentalData, plottingCondition,folderName);
    save parameters_Fit.mat x
    movefile('parameters_Fit.mat',folderName);
    
    %% Reporting real values removing Scaling
    realValue_Parameter(1)= x(1)/10;
    realValue_Parameter(2:4)=x(2:4)/1e3;
    realValue_Parameter
else
    %% Present results
%     x =[98.0000   23.0000   36.5000   53.5000];
%     x = [102.0000   23.0000   68.5000   53.5000];
    x= [106    22    50    66];

    plottingCondition=1;
    [experimentalData] = experimentalIntesityDistribution ;
    fitVal = fitData(x,experimentalData, plottingCondition,folderName);
    save parameters_Fit.mat x
    movefile('parameters_Fit.mat',folderName);
    
    % Reporting real values removing Scaling
    realValue_Parameter(1)= x(1)/10;
    realValue_Parameter(2:4)=x(2:4)/1e3;
    realValue_Parameter
end
cd ..