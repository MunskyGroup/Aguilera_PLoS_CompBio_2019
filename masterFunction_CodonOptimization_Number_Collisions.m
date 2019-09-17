clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces the mean and standard deviation for the number of collisions given in figures S4 to S6 in the publication.
% All results are displayed in command window.

%% Deffing parameter values
nRepetitions = 500;
nonConsideredInitialSimulationTime = 1000;
totalSimulationTime = 2000+ nonConsideredInitialSimulationTime;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;

%% loading experimental data and gene sequences
fileName = 'KDM5B.xls'; samplingRateInSeconds= 10; samplingForPlotting=1; maxlag =1000;geneFile = 'KDM5B_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.022;
% fileName = 'bact.xls'; samplingRateInSeconds= 3; samplingForPlotting=1; maxlag =300;  geneFile = 'Bactin_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.05;
% fileName = 'H2B.xls'; samplingRateInSeconds= 1; samplingForPlotting=1; maxlag =100;  geneFile = 'H2B_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.066;

%% Sequence optimizer
% optimization options 'optimized', 'deoptimized', 'natural', 'constant'
optimization = {'natural','optimized', 'deoptimized'};
for op =1: 3
    %% Creating Model from Gene sequence
    [generated_geneSequence, typeOfTag,tagPositions,parametersModel ,geneLength] = sequenceAnalyzer_Optimized(geneFile,k_initiation,k_elongationMean,optimization{op});
    %% Stochastic simulations
    parfor k = 1 : nRepetitions
             [~,~,collisions{k},~] = solve_direct_optimized_collisions(totalSimulationTime, parametersModel,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP)
             mean_collision(k) = mean(collisions{k});
    end
        mean_num_Collision(op) = mean(mean_collision);
        std_num_Collision(op) = std(mean_collision);
 
end

mean_num_Collision
std_num_Collision

cd ..
