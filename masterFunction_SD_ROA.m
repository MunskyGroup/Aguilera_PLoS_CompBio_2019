clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces figure S11 for the ROA experiments
%% All results are stored in Results_SD.

folderName = horzcat('Results_SD'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end

%% Reading the fasta format.
nonConsideredInitialSimulationTime =2000;
totalSimulationTime = 400+ nonConsideredInitialSimulationTime;
numberOfGenes=3;
timePerturbationApplication= nonConsideredInitialSimulationTime+200;
evaluatingInhibitor = 1;
evaluatingFRAP =0;
nPoints= totalSimulationTime;
nPoints2 = nPoints-nonConsideredInitialSimulationTime;
mean_intensityVector = zeros(numberOfGenes,nPoints2);
error_intensityVector= zeros(numberOfGenes,nPoints2);
selectedGene=2;
vector_NoRepetitions = [10,50,100];
meanHT_Delay=10;

for j=1:length(vector_NoRepetitions)
    nRepetitions = vector_NoRepetitions(j);
    s_out_Short {1,nRepetitions}=[];
    %for g =1:numberOfGenes
    for g =selectedGene:selectedGene
        
        if g == 1; geneFile = 'H2B_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.066; end
        if g == 2; geneFile = 'Bactin_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.05; end
        if g == 3; geneFile = 'KDM5B_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.022; end
        %% Creating Model from Gene sequence
        [generated_geneSequence, typeOfTag,tagPositions,parametersModel ,geneLength] = sequenceAnalyzer(geneFile,k_initiation,k_elongationMean);
        geneLength=geneLength-1;
        %% Stochastic simulations
        parfor k = 1 : nRepetitions
            timePerturbationApplication_withDelay = timePerturbationApplication + randn*meanHT_Delay;
            [t_out,s_out,~] = solve_direct_optimized(totalSimulationTime, parametersModel,timePerturbationApplication_withDelay,evaluatingInhibitor,evaluatingFRAP)
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
            preintensityVector = I;
            intensityVector (k,:)=preintensityVector;
        end
        
        meanVal =  mean (intensityVector(:,1));
        for k = 1 : nRepetitions
            intensityVector (k,:)= intensityVector (k,:)./meanVal;
        end
        
        
        time = linspace(0,totalSimulationTime-nonConsideredInitialSimulationTime,totalSimulationTime-nonConsideredInitialSimulationTime);
        
        %% calculating mean intensities and errors
        %         mean_intensityVector(g,:) = mean(intensityVector)./max(mean(intensityVector));
        mean_intensityVector(g,:) = mean(intensityVector);
        error_intensityVector(g,:) = std(intensityVector);%/sqrt(nRepetitions);
        meanSD=  mean (error_intensityVector(g,[1:200]) );
    end
    time=time-200;
    %% PLOTTING 1. Time Course
    figure('visible', 'off');
    fig1= gcf;
    fig1.PaperUnits = 'inches';
    fig1.PaperPosition = [0, 0, 2.2, 1.2];
    hold on
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
    text(50,1, ['SD =', num2str(round(meanSD,2))],'Color','red','FontSize',6);
    plot(zeros(1,10),linspace(0,1.5,10),'-r','LineWidth',1);
    xlim([0 totalSimulationTime]);
    xlabel('time (sec)','FontSize',14);
    ylabel('I(t) (a.u.)','FontSize',14);
    ylim([-0.05 1.5])
    xlim([-201 200]);
    box on
    set(gca,'linewidth',1)
    set (gca ,'FontSize',10, 'FontName', 'Arial');
    nameplot =[ 'HT_',num2str(nRepetitions)];
    print('-dpng','-r300',nameplot)
    movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
end

cd ..