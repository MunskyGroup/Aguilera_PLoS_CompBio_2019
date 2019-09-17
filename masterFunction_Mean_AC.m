clear all; close all; clc;
cd Source_Code
%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces figure 1 right bottom panel.
%% All results are stored in Results_Mean.
folderName = horzcat('Results_Mean'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end

%% Deffing parameter values
nRepetitions = 100;
nonConsideredInitialSimulationTime = 2000;
totalSimulationTime = 4000+ nonConsideredInitialSimulationTime;
k_initiation = 0.033;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;
selectedGene=3;

for g =1:1
    
    %% loading experimental data
    if g==1
        fileName = 'KDM5B.xls'; samplingRateInSeconds= 10; samplingForPlotting=1; maxlag =1000;geneFile = 'KDM5B_withTags.txt';k_elongationMean = 10.6; k_initiation=0.022;   fileName2 = 'KDM5B_AuCo.xls';
    elseif g==2
        fileName = 'bact.xls'; samplingRateInSeconds= 3; samplingForPlotting=1; maxlag =300;  geneFile = 'Bactin_withTags.txt'; k_elongationMean = 10.6;  k_initiation=0.05; fileName2 = 'ActB_AuCo.xls';
    elseif g==3
        fileName = 'H2B.xls'; samplingRateInSeconds= 1; samplingForPlotting=1; maxlag =100;  geneFile = 'H2B_withTags.txt'; k_elongationMean = 10.6; k_initiation=0.066; fileName2 = 'H2B_AuCo.xls';
    end
    
    %% Creating Model from Gene sequence
    [generated_geneSequence, typeOfTag,tagPositions,parametersModel ,geneLength] = sequenceAnalyzer(geneFile,k_initiation,k_elongationMean);
    geneLength=geneLength-1;
    %% Stochastic simulations
    parfor k = 1 : nRepetitions
        [t_out,s_out,RibosomePositions] = solve_direct_optimized(totalSimulationTime, parametersModel,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP);
        s_out_Short {1,k}= s_out(nonConsideredInitialSimulationTime+1:totalSimulationTime,:);
    end
    
    %% Creating intensity vectors
    intensityVector = zeros (nRepetitions,totalSimulationTime-nonConsideredInitialSimulationTime );
    I  = zeros (1,totalSimulationTime-nonConsideredInitialSimulationTime );
    probePositionVector = zeros (1, geneLength);
    for i =1: length(tagPositions)
        probePositionVector(1, tagPositions(i):end) = 1*i;
    end
    for k = 1 : nRepetitions
        for tp =1: totalSimulationTime-nonConsideredInitialSimulationTime
            temp_simulationOutPut = s_out_Short{1,k}(tp,1:geneLength);
            I(tp)= sum (probePositionVector .* temp_simulationOutPut );
        end
        intensityVector (k,:)= I;
    end
    intensityVector_nn= intensityVector; % not normalized vector
    
    %% time vector
    time = linspace(0,totalSimulationTime-nonConsideredInitialSimulationTime,totalSimulationTime-nonConsideredInitialSimulationTime);
    %% Calculating Autocorrelations for simulations.
    [lags, simulation_autocorrelation(g,:),simulation_sem_autocorrelation(g,:)] = AutoCorrelation_Simple_2(intensityVector,nRepetitions);
    %plotTimeCourses (time, intensityVector_nn,nRepetitions,geneFile,folderName)
end

%% Linear Fit
linreg = 1-((tagPositions(end)/2)/geneLength); %beginning of linear region
aT1 = find(simulation_autocorrelation(g,1:end)<linreg,1,'first');
aT2 = find(simulation_autocorrelation(g,1:end)<0.1,1,'first');
t_org =time;
[p,S] = polyfit(t_org([aT1:aT2]), simulation_autocorrelation(g,[aT1:aT2]),1);
[y_fit,~] = polyval(p,t_org([aT1:end]),S);
recovery_time= t_org(y_fit<0.05)+aT1;
recovery_time = recovery_time(1);
simulation_sem_autocorrelation = simulation_sem_autocorrelation.*0;


%% Plotting
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 2.2, 1];
lineProps.col = {[0.4 .0 1]};
lineProps.width = 2;
C = mseb(lags, simulation_autocorrelation(1,:),simulation_sem_autocorrelation(1,:),lineProps,0);
hold on;
plot ([recovery_time,recovery_time],[-0.05,1],'-k','LineWidth',0.5)
str = ['\tau_{FCS}'];
text(recovery_time+5,0.5,str,'Color','k','FontSize',8)
plot([-100,0,300],[0.05,0.05,0.05],':k','LineWidth',1) % ts
grid off;
box on
set(gca,'linewidth',1)
xlabel('\tau (sec)','FontSize',10);
ylabel('G(\tau)/G(0)','FontSize',10);
xlim([0 300])
%ylim([-0.1 1.4])
set (gca ,'FontSize',8, 'FontName', 'Arial');
nameplot ='AC';
print('-dpng','-r300',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
close all;




% figure('visible', 'off');
% fig1= gcf;
% fig1.PaperUnits = 'inches';
% fig1.PaperPosition = [0, 0, 2.2, 1];
% lineProps.col= {[1 0.6 0]};
% lineProps.width = 2;
% A = mseb(lags, simulation_autocorrelation(3,:),simulation_sem_autocorrelation(3,:),lineProps,0);
% lineProps.col={[0 .6 1]};
% lineProps.width = 2;
% B = mseb(lags, simulation_autocorrelation(2,:),simulation_sem_autocorrelation(2,:),lineProps,0);
% lineProps.col = {[0.4 .0 1]};
% lineProps.width = 2;
% C = mseb(lags, simulation_autocorrelation(1,:),simulation_sem_autocorrelation(1,:),lineProps,0);
% % hold on;
% % h=errorbar(experimental_lags(3,:), experimental_autocorrelation(3,:),experimental_sem(3,:), 'ok', 'MarkerEdgeColor',[0.1,0.1,0.1],'MarkerFaceColor',[1 0.6 0],'MarkerSize',3,'LineWidth',0.2);
% % h.CapSize = 0.1;
% % h=errorbar(experimental_lags(2,:), experimental_autocorrelation(2,:),experimental_sem(2,:), 'ok', 'MarkerEdgeColor',[0.1,0.1,0.1],'MarkerFaceColor',[0 .6 1],'MarkerSize',3,'LineWidth',0.2);
% % h.CapSize = 0.1;
% % h=errorbar(experimental_lags(1,:), experimental_autocorrelation(1,:),experimental_sem(1,:), 'ok', 'MarkerEdgeColor',[0.1,0.1,0.1],'MarkerFaceColor',[0.4 .0 1],'MarkerSize',3,'LineWidth',0.2);
% % h.CapSize = 0.1;
% lgd=legend([A.mainLine,B.mainLine,C.mainLine],'H2B','\beta-actin','KDM5B');
% set(lgd,'FontSize',7);
% lgd.Location='northeast';
% grid off;
% box on
% set(gca,'linewidth',1)
% xlabel('\tau (sec)','FontSize',10);
% ylabel('G(\tau)/G(0)','FontSize',10);
% xlim([0 300])
% %ylim([-0.1 1.4])
% set (gca ,'FontSize',8, 'FontName', 'Arial');
% nameplot ='AC';
% %print('-dpng','-r300',nameplot)
% print('-dpdf','-r300',nameplot)
% movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
% close all;

cd ..