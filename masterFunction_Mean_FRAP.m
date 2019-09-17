clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces figure 1 left bottom panel.
%% All results are stored in Results_Mean.

folderName = horzcat('Results_Mean'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end

%% Reading the fasta format.
nRepetitions = 100;
nonConsideredInitialSimulationTime =3000;
totalSimulationTime = 500+ nonConsideredInitialSimulationTime;
s_out_Short {1,nRepetitions}=[];
numberOfGenes=3;
timePerturbationApplication= nonConsideredInitialSimulationTime+100;

tp_r= 100;
evaluatingInhibitor = 0;
evaluatingFRAP =1;
nPoints= totalSimulationTime;
nPoints2 = nPoints-nonConsideredInitialSimulationTime;
mean_intensityVector = zeros(numberOfGenes,nPoints2);
error_intensityVector= zeros(numberOfGenes,nPoints2);

% for g =1:numberOfGenes
for g=3:3
    
    if g == 1; geneFile = 'H2B_withTags.txt';k_elongationMean = 10.6; k_initiation=0.066;  end
    if g == 2; geneFile = 'Bactin_withTags.txt'; k_elongationMean = 10.6;  k_initiation=0.05; end
    if g == 3; geneFile = 'KDM5B_withTags.txt';k_elongationMean = 10.6; k_initiation=0.022;end
    %% Creating Model from Gene sequence
    [generated_geneSequence, typeOfTag,tagPositions,parametersModel ,geneLength] = sequenceAnalyzer(geneFile,k_initiation,k_elongationMean);
    geneLength=geneLength-1;
    %% Stochastic simulations
    parfor k = 1 : nRepetitions
        [t_out,s_out,~,s_out_neg] = solve_direct_optimized(totalSimulationTime, parametersModel,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP)
        s_out_Short {1,k}= s_out(nonConsideredInitialSimulationTime+1:totalSimulationTime,:);
        s_out_Short_neg {1,k}= s_out_neg(nonConsideredInitialSimulationTime+1:totalSimulationTime,:);
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
            temp_simulationOutPut_neg = s_out_Short_neg{1,k}(tp,1:geneLength);
            I(tp)= sum (probePositionVector .* temp_simulationOutPut ) -...
                sum (probePositionVector .* temp_simulationOutPut_neg );
        end
        preintensityVector = I;
        intensityVector (k,:)=preintensityVector;
    end
    
    meanVal =  mean (intensityVector(:,end));
    for k = 1 : nRepetitions
        intensityVector (k,:)= intensityVector (k,:)./meanVal;
    end
    time = linspace(0,totalSimulationTime-nonConsideredInitialSimulationTime,totalSimulationTime-nonConsideredInitialSimulationTime);
    
    %% calculating mean intensities and errors
    mean_intensityVector(g,:) = mean(intensityVector);
    error_intensityVector(g,:) = std(intensityVector)./sqrt(nRepetitions);
    
end

t_org =time;
time=time-100;

%% Linear Fit
linreg = ((tagPositions(end)/2)/geneLength); %beginning of lionear region
aT1 = find(mean_intensityVector(g,tp_r+1:end)>linreg,1,'first');
aT2 = find(mean_intensityVector(g,tp_r+1:end)>0.9,1,'first');

[p,S] = polyfit(t_org([tp_r+aT1:1:tp_r+aT2]), mean_intensityVector(g,[tp_r+aT1:1:tp_r+aT2]),1);
[y_fit,~] = polyval(p,t_org([tp_r+aT1:1:end]),S);
recovery_time= t_org(y_fit>0.95)+aT1;
recovery_time = recovery_time(1);
ke  = (geneLength*0.95)./recovery_time;

%% PLOTTING 1. Time Course
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 2.2, 1];
hold on
% for g =1:3
for g=3:3
    if g == 1
        lineProps.col= {[1 0.6 0]};
        lineProps.width = 1.5;
        A = mseb(time,mean_intensityVector(g,:),error_intensityVector(g,:),lineProps,0); end
    if g == 2
        lineProps.col={[0 .6 1]};
        lineProps.width = 1.5;
        B = mseb(time,mean_intensityVector(g,:),error_intensityVector(g,:),lineProps,0); end
    if g == 3
        lineProps.col = {[0.4 .0 1]};
        lineProps.width = 1.5;
        C = mseb(time,mean_intensityVector(g,:),error_intensityVector(g,:),lineProps,0); end
    hold on
end

%% Plotting FRAP
plot(zeros(1,10),linspace(0,1.3,10),'-r','LineWidth',0.5);
plot([-100,0,500],[1,1,1],':k','LineWidth',1)
plot(t_org([tp_r+aT1:1:end])-100,y_fit,'-g','LineWidth',0.5)
plot(t_org(tp_r+aT1)-100,y_fit(1),'o','Color',[1, 0.8,0.2], 'MarkerSize',4,'MarkerFaceColor',[1, 0.8,0.2],'MarkerEdgeColor',[0.4 .0 1])
str1 = ['np_{FRAP}'];
text(t_org(tp_r+aT1)-100+22,0.1,str1,'Color','k','FontSize',6)
intersection = find(y_fit>0.95,1);
plot (recovery_time, y_fit(intersection),'o','Color',[1, 0.8,0.2], 'MarkerSize',4,'MarkerFaceColor',[1, 0.8,0.2],'MarkerEdgeColor',[0.4 .0 1])
plot ([recovery_time,recovery_time],[-0.05,1.3],':k','LineWidth',0.5)
str = ['\tau_{FRAP}'];
text(recovery_time+5,1.25,str,'Color','k','FontSize',6)
box on
set(gca,'linewidth',1)
xlim([-101 300]);
ylim([-0.05 1.4])
%grid on;
xlabel('time (sec)','FontSize',10);
ylabel('I(t) (au)','FontSize',10);
set (gca ,'FontSize',8, 'FontName', 'Arial');
nameplot = 'FRAP';
print('-dpng','-r300',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');

cd ..
