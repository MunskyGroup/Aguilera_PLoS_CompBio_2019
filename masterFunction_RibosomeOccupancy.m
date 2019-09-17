clear all; close all;clc;
cd Source_Code

%% This code is inteded to display the sequence composition for the studied genes.
%% The code reproduces figure S3 .
%% All results are stored in Results_CAI.

folderName = horzcat('Results_CAI'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end

%% Deffing parameter values
nRepetitions = 100;
nonConsideredInitialSimulationTime = 10000;
totalSimulationTime = 5000+ nonConsideredInitialSimulationTime;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;

for g =1:3
%% loading experimental data
    if g == 1; geneFile = 'H2B_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.066;end
    if g == 2; geneFile = 'Bactin_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.05;end
    if g == 3; geneFile = 'KDM5B_withTags.txt'; k_elongationMean = 10.6;k_initiation = 0.022;end
%% Creating Model from Gene sequence
[generated_geneSequence, typeOfTag,tagPositions,parametersModel ,geneLength] = sequenceAnalyzer(geneFile,k_initiation,k_elongationMean);
lengthTagRegion=tagPositions(end);
geneLength=geneLength-1;
%% Stochastic simulations
 s_out {nRepetitions}=[];
pre_ribosomeDensity = zeros(nRepetitions, geneLength);
parfor k = 1 : nRepetitions
    [~,s_out{k},~,~] = solve_direct_optimized(totalSimulationTime, parametersModel, timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP);
 pre_ribosomeDensity (k,:) = sum (s_out{k}(nonConsideredInitialSimulationTime+1:totalSimulationTime,:));
end
ribosomeDensity = mean(pre_ribosomeDensity);
ribosomeDensity = smooth(ribosomeDensity,9);
sum_ribosomeDensity  = sum(ribosomeDensity);
norm_ribosomeDensity = ribosomeDensity./sum_ribosomeDensity;
bins = linspace(1,geneLength,geneLength);
%% PLOTTING 
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 2.3, 1.2];
if  strcmp(geneFile , 'H2B_withTags.txt') == 1
bar(bins,norm_ribosomeDensity,'FaceColor',[1 .6 .0],'EdgeColor',[1 .6 .0])
hold on
bar(bins(1:lengthTagRegion),norm_ribosomeDensity(1:lengthTagRegion),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
end
if strcmp(geneFile ,'Bactin_withTags.txt') == 1
bar(bins,norm_ribosomeDensity,'FaceColor',[0 .6 1],'EdgeColor',[0 .6 1])
hold on
bar(bins(1:lengthTagRegion),norm_ribosomeDensity(1:lengthTagRegion),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
end
if strcmp(geneFile ,'KDM5B_withTags.txt') == 1 
bar(bins,norm_ribosomeDensity,'FaceColor',[0.4 .0 1],'EdgeColor',[0.4 .0 1])
hold on
bar(bins(1:lengthTagRegion),norm_ribosomeDensity(1:lengthTagRegion),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
end
%xlim([1 1900])
%ylim([0 1])
box on
set(gca,'linewidth',1)
xlabel('Codon Position','FontSize',12);
ylabel('P(Ribosome)','FontSize',12);
set (gca ,'FontSize',8, 'FontName', 'Arial');
nameplot = horzcat('RO_',geneFile(1:end-13));
print('-dpng','-r600',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
end
cd ..