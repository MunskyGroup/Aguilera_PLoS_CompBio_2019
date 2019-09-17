clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces figure S11 for the FRAP experiments
%% All results are stored in Results_SD.

folderName = horzcat('Results_SD'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end
%% Reading the fasta format.
nonConsideredInitialSimulationTime =3000;
totalSimulationTime = 400+ nonConsideredInitialSimulationTime;
numberOfGenes=3;
timePerturbationApplication= nonConsideredInitialSimulationTime+200;
tp_r= 200;
evaluatingInhibitor = 0;
evaluatingFRAP =1;
nPoints= totalSimulationTime;
nPoints2 = nPoints-nonConsideredInitialSimulationTime;
mean_intensityVector = zeros(numberOfGenes,nPoints2);
error_intensityVector= zeros(numberOfGenes,nPoints2);
selectedGene=2;
vector_NoRepetitions = [10,50,100];
for j=1:length(vector_NoRepetitions)
    nRepetitions = vector_NoRepetitions(j);
    s_out_Short {1,nRepetitions}=[];
    %for g =1:numberOfGenes
    for g=selectedGene:selectedGene
        if g == 1; geneFile = 'H2B_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.066; end
        if g == 2; geneFile = 'Bactin_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.05; end
        if g == 3; geneFile = 'KDM5B_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.022;end
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
            intensityVector (k,:)=I;
        end
        
        meanVal =  mean (intensityVector(:,end));
        for k = 1 : nRepetitions
            intensityVector (k,:)= intensityVector (k,:)./meanVal;
        end
        
        time = linspace(0,totalSimulationTime-nonConsideredInitialSimulationTime,totalSimulationTime-nonConsideredInitialSimulationTime);
        %% calculating mean intensities and errors
        mean_intensityVector(g,:) = mean(intensityVector);
        error_intensityVector(g,:) = std(intensityVector);;
        meanSD= mean (error_intensityVector(g,[end-100:end]));
    end
    time=time-200;
    %% PLOTTING 1. Time Course
    figure('visible', 'off');
    fig1= gcf;
    fig1.PaperUnits = 'inches';
    fig1.PaperPosition = [0, 0, 2.2, 1.2];
    hold on
    %for g =1:3
    for g=selectedGene:selectedGene
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
    text(10,2.5, ['SD =', num2str(round(meanSD,2))],'Color','red','FontSize',6);
    plot(zeros(1,10),linspace(0,3,10),'-r','LineWidth',1);
    box on
    set(gca,'linewidth',1)
    xlim([-201 200]);
    ylim ([0,3])
    xlabel('time (sec)','FontSize',14);
    ylabel('I(t) (a.u.)','FontSize',14);
    set (gca ,'FontSize',10, 'FontName', 'Arial');
    nameplot =[ 'FRAP_',num2str(nRepetitions)];
    print('-dpng','-r300',nameplot)
    movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
end
cd ..