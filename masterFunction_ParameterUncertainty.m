close all; clc; clear all;
cd Source_Code

folderName = horzcat('Results_ParameterUncertainty_final'); if exist (folderName, 'dir') ~= 7; mkdir(folderName); end

%% Code to calculate parameter uncertanty
%% This code is inteded to simulate single molecule stochastic dynamics and calculate parameter uncertanty
%% The code reproduces the SD for the parameters given in the table of parameters.
%% All results are directly printed to screen.

%% Deffining the previously optmized parameters
% To facilitate the uncertainty process the parameters are scaled to
% integer values. Elongation rates are scaled as ke_scale = ke_real *10
% Initiation rates are scaled as ki_scaled = ki_real*1e4;
pars_best =  [106    22    50    66];
pars_real = pars_best;
numberTotalEvaluations=1000;

%ke  = pars_new(1); % Elongation rate
%ki_gene1 = pars_new(2);  %% Intitiation rate for KDM5B
%ki_gene2 =pars_new(3);   %% Intitiation rate for B-act
%ki_gene3 = pars_new(4);  %% Intitiation rate for H2B

% load experimental data
[renormalized_intensities] = experimentalIntesityDistribution ;

%% Using fast SSA. Hybrid stochastic before the FSS and deterministic after FSS.
finvalmin=inf;
ftvals_best = [10 10 10 10];
pars = pars_best;
pars_sv = [];
inorout = [];
J_sv = [];
Ftvals = [1 1 1 1 ];
for ip=1:numberTotalEvaluations
    ip
    if ip==1
        delt = 0;
    else
        delt = 0.1;
    end
    %% Deffing parameters.
    pars_new = pars.*(1+delt*randn(size(pars)));
    pars_sv = [pars_sv;pars_new];
    finval = fitData(pars_new,renormalized_intensities, 0,'');
    % Updates best parametes if it finds something better than given best parameters
    if finval<finvalmin ==1
        pars_best = pars_new;
        ftvals_best =  finval;
        finvalmin = finval;
    end
    
    if finval <= 1.1*ftvals_best ==1
        pars = pars_new;
        inorout = [inorout,1];
        J_sv = [J_sv;Ftvals];
    else
        inorout = [inorout,0];
    end
    
end
Parameter_AboveThreshold = pars_sv(inorout==1,:);

%% Save best values
cd (folderName)
S = dir('randomSearchData*.*');
try
    numberOf_RS_Repetitions =  length (S);
catch
    numberOf_RS_Repetitions = 0;
end
fname=[pwd,'/randomSearchData_',num2str(numberOf_RS_Repetitions+1),'.mat'];
save(fname,'Parameter_AboveThreshold','pars_best','J_sv')

%% Load previous list of parameters that fullfil the selection criterion
sel_param_list = [0,0,0,0];
% load all repetitions to Plot data.
S_load = dir('randomSearchData*.*');
numberOf_RS_Repetitions =  length (S_load);
for i = 1: numberOf_RS_Repetitions
    fileNames =  ['randomSearchData_',num2str(i),'.mat'];
    load (fileNames)
    sel_param_list= [sel_param_list;Parameter_AboveThreshold];
end
% removing rows with only zeros.
sel_param_list = sel_param_list(any(sel_param_list,2),:);

realValue_Parameter_AboveThreshold(:,1)= Parameter_AboveThreshold(:,1)/10;
realValue_Parameter_AboveThreshold(:,2:4)=Parameter_AboveThreshold(:,2:4)/1e3;
mean (realValue_Parameter_AboveThreshold)
std (realValue_Parameter_AboveThreshold)

cd ..

