function plotCombinationDistribution(geneFile,meanIntensity_simulation,folderName)

%% Plotting Distribution for HA sequence
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 1.5, 1.];
hold on
if  strcmp(geneFile , 'H2B_withTags.txt') == 1
    colorPlot =[1 0.6 0];
    g=3;
end
if strcmp(geneFile ,'Bactin_withTags.txt') == 1
    colorPlot =[0 .6 1];
    g=2;
end
if strcmp(geneFile ,'KDM5B_withTags.txt') == 1
    colorPlot =[0.4 .0 1];
    g=1;
end
lineProps.width = 1;
edg= 0:1:12;

%% this section calculates the Distributions for the experimental data

%% For the error bar use 300 spots.
%% For the solid line use total data.

nRep=100;
for i=1: nRep
    [randIntensities,AllIntensity] = randomSampleOfCells(g,nRep);
    hisData_sampling(i,:) = histcounts (randIntensities{i},edg,'Normalization','probability');
    [xb,Data_stairs(:,i)] = stairs(hisData_sampling(i,:));
end
Data_stairs= Data_stairs';
std_stairs_ExpData =  std(Data_stairs);



hisData_All = histcounts (AllIntensity,edg,'Normalization','probability');
[xb,Data_stairs_All] = stairs(hisData_All);
All_stairs_expData = Data_stairs_All';

lineProps.col= {[0,0,0]};
h = mseb(xb', All_stairs_expData,std_stairs_ExpData,lineProps,0);


%% this section calculates the Distributions for the simulateions 
sampleSize =max(size(meanIntensity_simulation));
for i=1: nRep
        ii = randi([1, sampleSize], 1, sampleSize);
        sample_meanIntensity_simulation = meanIntensity_simulation(ii);
        hisSim_sampling(i,:) = histcounts (sample_meanIntensity_simulation,edg,'Normalization','probability');
      [xb,sim_stairs(:,i)] = stairs(hisSim_sampling(i,:));
end
sim_stairs= sim_stairs';
mean_stairs_Sim = mean(sim_stairs);
std_stairs_Sim =  std(sim_stairs);
lineProps.col= {colorPlot};
h_sim = mseb(xb', mean_stairs_Sim,std_stairs_Sim,lineProps,0);
xlim([1 , 14])
ylim([0 ,0.6])
lgd =legend([h.mainLine, h_sim.mainLine],'Data','Model');
% lgd =legend([h, h_sim],'Data','Model');
set(lgd,'FontSize',4);
lgd.Location='northeast';

box on
set(gca,'linewidth',1)
grid off
set (gca ,'FontSize',6);
set(gca, 'FontName', 'Arial');
xlabel('Intensity (ump)','FontName', 'Arial','FontSize',8);
ylabel('Probability','FontName', 'Arial','FontSize',8);
nameplot = horzcat('Dist_Comp_',geneFile(1:end-13));
print('-dpng','-r300',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
end