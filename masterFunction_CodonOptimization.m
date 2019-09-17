clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces figures S4 to S6 in the publication.
% All results are stored in Results_Codon_Optimization.

folderName = horzcat('Results_Codon_Optimization'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end

%% Deffing parameter values
nRepetitions = 500;
nonConsideredInitialSimulationTime = 1000;
totalSimulationTime = 2000+ nonConsideredInitialSimulationTime;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;

%% loading experimental data and gene sequences
fileName = 'KDM5B.xls'; samplingRateInSeconds= 10; samplingForPlotting=1; maxlag =1000;geneFile = 'KDM5B_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.022;
% fileName = 'bact.xls'; samplingRateInSeconds= 3; samplingForPlotting=1; maxlag =300;  geneFile = 'Bactin_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.05;
% fileName = 'H2B.xls'; samplingRateInSeconds= 1; samplingForPlotting=1; maxlag =100;  geneFile = 'H2B_withTags.txt'; k_elongationMean = 10.6; k_initiation = 0.066;

%% Sequence optimizer
% optimization options 'optimized', 'deoptimized', 'natural', 'constant'
optimization = {'natural','optimized', 'deoptimized'};
for op =1: 3
    %% Creating Model from Gene sequence
    [generated_geneSequence, typeOfTag,tagPositions,parametersModel ,geneLength] = sequenceAnalyzer_Optimized(geneFile,k_initiation,k_elongationMean,optimization{op});
    %% Stochastic simulations
    parfor k = 1 : nRepetitions
        [t_out,s_out{k},~,~] = solve_direct_optimized(totalSimulationTime, parametersModel,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP);
        s_out_Short {1,k}= s_out{k}(nonConsideredInitialSimulationTime+1:totalSimulationTime,:);
    end
    
    for k = 1 : nRepetitions
        pre_ribosomeDensity (k,:) = sum (s_out{k}(nonConsideredInitialSimulationTime+1:totalSimulationTime,:));
    end
    
    %% Ribosome loading
    ribosomeDensity = mean(pre_ribosomeDensity);
    ribosomeDensity = smooth(ribosomeDensity,9);
    sum_ribosomeDensity  = sum(ribosomeDensity);
    norm_ribosomeDensity{op} = ribosomeDensity./sum_ribosomeDensity;
    bins_RibosomeLoading{op} = linspace(1,geneLength,geneLength);
    
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
    
    %% Intensity Distributions
    hisData (op,:)=  intensityVector(:,end);
    %% Calculate autocorrelation
    [lags, simulation_autocorrelation(op,:),~] = AutoCorrelation_Simple_2(intensityVector,nRepetitions);
end
time = linspace(0,totalSimulationTime-nonConsideredInitialSimulationTime,totalSimulationTime-nonConsideredInitialSimulationTime);

%% Plotting histograms
plottingIntensityDistribution(hisData,geneFile, folderName);

%% PLOTTING  Ribosome loading
colorPlot_OP= [1 0 0];
colorPlot_DO= [0 0 0];
for i =1:3
    figure('visible', 'off');
    fig1= gcf;
    fig1.PaperUnits = 'inches';
    fig1.PaperPosition = [0, 0, 2.4, 0.6];
    if  strcmp(geneFile , 'H2B_withTags.txt') == 1
        switch i
            case 1
                plot (bins_RibosomeLoading{i},norm_ribosomeDensity{i},'-','Color',[1 0.6 0], 'LineWidth',1 )
            case 2
                plot (bins_RibosomeLoading{i},norm_ribosomeDensity{i},'-','Color',colorPlot_OP, 'LineWidth',1 )
            case 3
                plot (bins_RibosomeLoading{i},norm_ribosomeDensity{i},'-','Color',colorPlot_DO, 'LineWidth',1 )
        end
        xticks([0 200 400 ])
        xticklabels({'0','200','400'})
        
    end
    if strcmp(geneFile ,'Bactin_withTags.txt') == 1
        switch i
            case 1
                plot (bins_RibosomeLoading{i},norm_ribosomeDensity{i},'-', 'Color',[0 .6 1], 'LineWidth',1 )
            case 2
                plot (bins_RibosomeLoading{i},norm_ribosomeDensity{i},'-', 'Color',colorPlot_OP, 'LineWidth',1 )
            case 3
                plot (bins_RibosomeLoading{i},norm_ribosomeDensity{i},'-', 'Color',colorPlot_DO, 'LineWidth',1 )
        end
        xticks([0 200 400 600])
        xticklabels({'0','200','400','600'})
    end
    if strcmp(geneFile ,'KDM5B_withTags.txt') == 1
        switch i
            case 1
                plot (bins_RibosomeLoading{i},norm_ribosomeDensity{i},'-','Color',[0.4 .0 1], 'LineWidth',1 )
            case 2
                plot (bins_RibosomeLoading{i},norm_ribosomeDensity{i},'-','Color',colorPlot_OP, 'LineWidth',1 )
            case 3
                plot (bins_RibosomeLoading{i},norm_ribosomeDensity{i},'-','Color',colorPlot_DO, 'LineWidth',1 )
        end
        xticks([0 500 1000 1500 ])
        xticklabels({'0','500','1000','1500'})
        
    end
    
    maxY= max([max(norm_ribosomeDensity{1}),max(norm_ribosomeDensity{2}),max(norm_ribosomeDensity{3})]);
    
    xlim([1 max(bins_RibosomeLoading{i})])
    ylim([0 maxY])
    box on
    set(gca,'linewidth',1)
    xlabel('Ribosome Position','FontSize',8);
    set (gca ,'FontSize',6, 'FontName', 'Arial');
    nameplot = horzcat('TC_',geneFile(1:end-13),'_',num2str(i) );
    print('-dpng','-r300',nameplot)
    movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
end

%% Plotting
colorPlot_OP= [1 0 0];
colorPlot_DO= [0 0 0];

for i =1:3
    simulation_autocorrelation(i,:) = simulation_autocorrelation(i,:)+ abs(min (simulation_autocorrelation(i,:)));
    simulation_autocorrelation(i,:) = simulation_autocorrelation(i,:)/max(simulation_autocorrelation(i,:));
end

space = 9;
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 2.5, 1.6];
hold on
if  strcmp(geneFile , 'H2B_withTags.txt') == 1
    a=plot (time(1:space:end),simulation_autocorrelation(1,1:space:end),'-', 'Color',[1 0.6 0], 'LineWidth',1);
    b=plot (time(1:space:end),simulation_autocorrelation(2,1:space:end),':*', 'Color',colorPlot_OP, 'LineWidth',1,'MarkerFaceColor',[1 0.6 0],'MarkerSize',4);
    c=plot (time(1:space:end),simulation_autocorrelation(3,1:space:end),'-.s', 'Color',colorPlot_DO, 'LineWidth',1,'MarkerFaceColor',[1 0.6 0],'MarkerSize',3);
    xlim([0 120])
    
end
if  strcmp(geneFile , 'Bactin_withTags.txt') == 1
    a=plot (time(1:space:end),simulation_autocorrelation(1,1:space:end),'-', 'Color',[0 .6 1], 'LineWidth',2);
    b=plot (time(1:space:end),simulation_autocorrelation(2,1:space:end),':*', 'Color',colorPlot_OP, 'LineWidth',2,'MarkerFaceColor',[0 .6 1],'MarkerSize',4);
    c=plot (time(1:space:end),simulation_autocorrelation(3,1:space:end),'-.s', 'Color',colorPlot_DO, 'LineWidth',2,'MarkerFaceColor',[0 .6 1],'MarkerSize',3);
    xlim([0 120])
    
end
if  strcmp(geneFile , 'KDM5B_withTags.txt') == 1
    a=plot (time(1:space:end),simulation_autocorrelation(1,1:space:end),'-', 'Color',[0.4 .0 1], 'LineWidth',1);
    b=plot (time(1:space:end),simulation_autocorrelation(2,1:space:end),':*', 'Color',colorPlot_OP, 'LineWidth',1,'MarkerFaceColor',[0.4 .0 1],'MarkerSize',4);
    c=plot (time(1:space:end),simulation_autocorrelation(3,1:space:end),'-.s', 'Color',colorPlot_DO, 'LineWidth',1,'MarkerFaceColor',[0.4 .0 1],'MarkerSize',3);
    xlim([0 350])
end
grid off;
box on
set(gca,'linewidth',1)
ylim([-0.1 1])
xlabel('\tau (sec)','FontSize',12);
ylabel('G(\tau)/G(0)','FontSize',12);
nameplot = horzcat('Opt_',geneFile(1:end-13));
lgd=legend([a,b,c],'Natural','Common','Rare');
set(lgd,'FontSize',8);
lgd.Location='northeast';
set (gca ,'FontSize',8, 'FontName', 'Arial');
print('-dpng','-r300',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
close all;

cd ..
