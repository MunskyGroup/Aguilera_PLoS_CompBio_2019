clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces figure 5 bottom panel in the publication.
%% All results are stored in Results_Collisions_Sweep
%% For the Paper, b-act gene is selected.

folderName = horzcat('Results_Collisions_Sweep'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end
%% Deffing parameter values
nRepetitions = 10;
nonConsideredInitialSimulationTime = 2000;
totalSimulationTime = 2000+ nonConsideredInitialSimulationTime;
Resolution =30;
geneLength =1000;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;
k_elongationMean_vector = linspace(2,13,Resolution);
k_initiation_vector = linspace(0.005,0.08,Resolution);
usingSynteticGene =1;

%% Running Simulations
for i =1:length(k_elongationMean_vector)
    for j =1:length(k_initiation_vector)
        k_elongationMean = k_elongationMean_vector(1,i);
        k_initiation = k_initiation_vector(1,j);
        if usingSynteticGene==1
            k_elongation = ones(1,geneLength-2)*k_elongationMean;
            parametersModel{j} = [k_initiation,k_elongation,10];
        else
            geneFile = 'Bactin_withTags.txt';
            [~, ~,~,parametersModel{j}] = sequenceAnalyzer(geneFile,k_initiation,k_elongationMean);
        end
        parfor k = 1 : nRepetitions
            [~,~,collisions{k},~] = solve_direct_optimized_collisions(totalSimulationTime, parametersModel{j},timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP)
            mean_collision(k) = mean(collisions{k});
            median_collision(k) = median(collisions{k});
        end
        probability_Collision_mean(i,j) = mean(mean_collision);
        probability_Collision_median(i,j) = mean(median_collision); % notice that the median is calculated inside the parfor.
    end
end

%% saving data
save('collision_data.mat' ,'collisions','probability_Collision_mean', 'probability_Collision_median','k_elongationMean_vector','k_initiation_vector')
%% Plotting
if usingSynteticGene ==1
    nameplot ='mean_Collision';
else
    nameplot =['mean_Collision_',geneFile([1:6])];
end
plot_Collisions (probability_Collision_mean,k_initiation_vector,k_elongationMean_vector,nameplot,folderName)

cd ..