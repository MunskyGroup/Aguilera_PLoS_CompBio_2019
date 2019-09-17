clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces figure S2 in the publication.
% All results are stored in Results_Collisions.

folderName = horzcat('Results_Collisions'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end
%% Deffing parameter values
nonConsideredInitialSimulationTime = 1000;
totalSimulationTime = 10000+ nonConsideredInitialSimulationTime;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;
k_elongationMean_vector = [2,4,10];
k_initiation_vector = [0.06];
usingSynteticGene =0;
geneLength =1000;

for i =1:length(k_initiation_vector)
for j =1:length( k_elongationMean_vector )
k_elongationMean = k_elongationMean_vector(1,j);
k_initiation = k_initiation_vector(1,i);
if usingSynteticGene==1
k_elongation = ones(1,geneLength-2)*k_elongationMean;
parametersModel{j} = [k_initiation,k_elongation,10];
else
geneFile = 'Bactin_withTags.txt'; 
[~, ~,~,parametersModel{j}] = sequenceAnalyzer(geneFile,k_initiation,k_elongationMean);
end
[~,s_out{j},collisions{j},RibosomeDynamcis{j}] = solve_direct_optimized_collisions(totalSimulationTime, parametersModel{j},timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP);
end
end

%% saving data
%save('collision_data_new.mat','s_out','collisions', 'k_elongationMean_vector','k_initiation_vector','geneLength')

%% Plotting Distribution of Collisions
for j =1:length(k_elongationMean_vector)
 nameplot =['ke_', strrep(['_',num2str(k_elongationMean_vector(j))],'.','_')]; 
 plotCollisionDistributions_ke (collisions{j},nameplot,folderName)
%  max (collisions{j})
end

%% Plotting Kymographs 
for j =1:length(k_elongationMean_vector)
 nameplot =strrep(['ke_kymograph_',num2str(k_elongationMean_vector(j))],'.','_');
 kymograph (s_out{1,j},nameplot,folderName);
end

cd ..