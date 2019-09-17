clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces figure S11 for the AC experiments
%% All results are stored in Results_SD.

folderName = horzcat('Results_SD'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end
%% Deffing parameter values
nonConsideredInitialSimulationTime = 2000;
totalSimulationTime = 2000+ nonConsideredInitialSimulationTime;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;
selectedGene=2;
vector_NoRepetitions = [10,50,100];

for j=1:length(vector_NoRepetitions)
    nRepetitions = vector_NoRepetitions(j);
    
    %for g =1:3
    for g =selectedGene:selectedGene
        
        %% loading experimental data
        if g==1
            fileName = 'KDM5B.xls'; samplingRateInSeconds= 10; samplingForPlotting=1; maxlag =1000;geneFile = 'KDM5B_withTags.txt'; k_elongationMean = 10.6;k_initiation = 0.022;
        elseif g==2
            fileName = 'bact.xls'; samplingRateInSeconds= 3; samplingForPlotting=1; maxlag =300;  geneFile = 'Bactin_withTags.txt'; k_elongationMean = 10.6;k_initiation = 0.05;
        elseif g==3
            fileName = 'H2B.xls'; samplingRateInSeconds= 1; samplingForPlotting=1; maxlag =100;  geneFile = 'H2B_withTags.txt'; k_elongationMean = 10.6;k_initiation = 0.066;
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
        %% time vector
        time = linspace(0,totalSimulationTime-nonConsideredInitialSimulationTime,totalSimulationTime-nonConsideredInitialSimulationTime);
        %% Calculating Autocorrelations for simulations.
        [lags, simulation_mean_autocorrelation(g,:),simulation_sem_autocorrelation(g,:)] = AutoCorrelation_Simple_2(intensityVector,nRepetitions);
        simulation_sd_autocorrelation (g,:) =  simulation_sem_autocorrelation(g,:)* sqrt(nRepetitions);
        meanSD= mean (simulation_sd_autocorrelation(g,:));
        
    end
    
    %% Plotting 1 - Autocorrelation.
    figure('visible', 'off');
    fig1= gcf;
    fig1.PaperUnits = 'inches';
    fig1.PaperPosition = [0, 0, 2.2, 1.2];
    switch selectedGene
        case 1
            lineProps.col= {[1 0.6 0]};
            lineProps.width = 2;
            A = mseb(lags, simulation_mean_autocorrelation(3,:),simulation_sd_autocorrelation(3,:),lineProps,0);
        case 2
            lineProps.col={[0 .6 1]};
            lineProps.width = 2;
            B = mseb(lags, simulation_mean_autocorrelation(2,:),simulation_sd_autocorrelation(2,:),lineProps,0);
        case 3
            lineProps.col = {[0.4 .0 1]};
            lineProps.width = 2;
            C = mseb(lags, simulation_mean_autocorrelation(1,:),simulation_sd_autocorrelation(1,:),lineProps,0);
    end
    
    text(150,0.8, ['SD =', num2str(round(meanSD,2))],'Color','red','FontSize',6);
    
    hold on;
    box on
    set(gca,'linewidth',1)
    xlabel('time delay \tau (sec)','FontSize',14);
    ylabel('G(\tau)/G(0)','FontSize',14);
    xlim([0 350])
    % ylim([-0.5 1.3])
    ylim([-0.5 1.3])
    set (gca ,'FontSize',10, 'FontName', 'Arial');
    nameplot =[ 'ac_SD_',num2str(nRepetitions)];
    print('-dpng','-r300',nameplot)
    movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
end
cd ..