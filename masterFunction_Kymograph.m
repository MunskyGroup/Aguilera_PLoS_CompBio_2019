clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces figure 2 left bottom panel.
%% All results are stored in Results_Kymograph.

folderName = horzcat('Results_Kymograph'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end
%% Deffing parameter values
nonConsideredInitialSimulationTime = 0;
totalSimulationTime = 1000+ nonConsideredInitialSimulationTime;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;
tagPositions = [2,11,20,196,206,218,228,300,309,318];
nRepetitions=1;
k_elongationMean_vector = 9.8;
k_initiation_vector = 0.039;
usingSynteticGene =1;
geneLength =1000;
geneLength = geneLength+tagPositions(end);

%% Running stochastic simulations
for i =1:length(k_elongationMean_vector)
    parfor j =1:length(k_initiation_vector)
        k_elongationMean = k_elongationMean_vector(1,i);
        k_initiation = k_initiation_vector(1,j);
        if usingSynteticGene==1
            k_elongation = ones(1,geneLength-2)*k_elongationMean;
            parametersModel{j} = [k_initiation,k_elongation,10];
        else
            geneFile = 'Bactin_withTags.txt';
            [~, ~,~,parametersModel{j} ,~] = sequenceAnalyzer(geneFile,k_initiation,k_elongationMean);
        end
        [~,s_out{j},collisions{j},RibosomeDynamcis{j}] = solve_direct_optimized_collisions(totalSimulationTime, parametersModel{j},timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP);
    end
end

%% Creating intensity vectors
intensityVector = zeros (1,totalSimulationTime-nonConsideredInitialSimulationTime );
I  = zeros (1,totalSimulationTime-nonConsideredInitialSimulationTime );
probePositionVector = zeros (1, geneLength-2);
for i =1: length(tagPositions)
    probePositionVector(1, tagPositions(i):end) = 1*i;
end
for k = 1 : nRepetitions
    for tp =1: totalSimulationTime-nonConsideredInitialSimulationTime
        temp_simulationOutPut = s_out{1}(tp,1:geneLength-2);
        I(tp)= sum (probePositionVector .* temp_simulationOutPut );
    end
    preintensityVector = I;
    intensityVector (k,:)= preintensityVector./max(preintensityVector);
end
pre_ribosomeDensity= sum (s_out{1}(:,:),2)';

%% Converting the simulation results to a binary matrix
s_out_kym =s_out{1}(1:1000,:);
for i =1:size(s_out_kym,1)
    for j =2:size(s_out_kym,2)-10
        if s_out_kym(i,j-1) ==0 && s_out_kym(i,j)==1
            s_out_kym(i,j:j+10)=1;
        end
    end
end

% Plotting Kymograph
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 4, 2.5];
hold on
imshow(s_out_kym)
axis tight
h = gca;
h.Visible = 'On';
box on
set(gca,'linewidth',2)
set (gca ,'FontSize',12);
set(gca, 'FontName', 'Arial');
%xlabel('Codon Position','FontSize',14,'FontName', 'Arial');
set(gca,'xtick',[])
set(gca,'yticklabel',[])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% ylabel('Time (sec)','FontName', 'Arial','FontSize',14);
nameplot = 'desc_kym';
print('-dpng','-r600',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');


%% plotting time vs intensity
time = linspace(0,1000,1000);
selectedSpot=1;
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0,  2.9, 0.5];
%plot(time,intensityVector(selectedSpot,:),'Color',[0.2 0.2 0.2],'LineWidth',1);
rectangle('Position',[1,0,1000,15],'FaceColor',[0.1 0.1 0.1])
hold on
% plot(time,pre_ribosomeDensity,'Color',[0.2 0.2 0.2],'LineWidth',1);
plot(time,pre_ribosomeDensity,'Color',[0.9 0.9 0.9],'LineWidth',1);
xlim([0 1000])
%ylim([-0.5 1])
box on
set(gca,'linewidth',1)
%xlabel('time (sec)','FontSize',12);
%ylabel('I(t) (a.u.)','FontSize',12);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set (gca ,'FontSize',8, 'FontName', 'Arial');
nameplot = 'desc_TC';
print('-dpng','-r600',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');

cd ..

