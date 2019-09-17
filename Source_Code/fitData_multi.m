%% This code is inteded to simulate single molecule stochastic dynamics.
function fitVal = fitData_multi(x,renormalized_intensities, plottingCondition,folderName)

x=round(x);



k_elongation1 = x(1)/10;
% k_elongation2 = x(2)/10;
% k_elongation3 = x(3)/10;

k_initiation_1 = x(2)/1000;
k_initiation_2 = x(3)/1000;
k_initiation_3 = x(4)/1000;


%% Deffing parameter values
nRepetitions = 100;
nonConsideredInitialSimulationTime = 1000;
totalSimulationTime = 2000+ nonConsideredInitialSimulationTime;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;
%% loading experimental data
c=4;
for g=1:3
    c=c-1;
    if g==1
        fileName = 'KDM5B.xls'; samplingRateInSeconds= 10; maxlag =250;geneFile = 'KDM5B_withTags.txt'; fileName_AC = 'KDM5B_AuCo.xls'; ExpRib=4.99; k_initiation =k_initiation_1; k_elongationMean=k_elongation1;
    end
    if g==2
        fileName = 'bact.xls'; samplingRateInSeconds= 3; maxlag =150;  geneFile = 'Bactin_withTags.txt';  fileName_AC = 'ActB_AuCo.xls'; ExpRib=3.01; k_initiation =k_initiation_2;k_elongationMean=k_elongation1;
    end
    if g==3
        fileName = 'H2B.xls'; samplingRateInSeconds= 1; maxlag =75;  geneFile = 'H2B_withTags.txt'; fileName_AC = 'H2B_AuCo.xls'; ExpRib=2.12; k_initiation=k_initiation_3;k_elongationMean=k_elongation1;
    end

    %% Creating Model from Gene sequence
    [~, ~,tagPositions,parametersModel ,geneLength] = sequenceAnalyzer(geneFile,k_initiation,k_elongationMean);
    geneLength=geneLength-1;
    %% Stochastic simulations
    parfor k = 1 : nRepetitions
        [~,s_out,~,~] = solve_direct_optimized(totalSimulationTime, parametersModel,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP);
        s_out_Short {1,k}= s_out(nonConsideredInitialSimulationTime+1:totalSimulationTime,:);
        pre_ribosomeDensity(k,:) = sum (s_out(nonConsideredInitialSimulationTime+1:totalSimulationTime,:),2)';
    end
    NoRibosomes = mean(pre_ribosomeDensity(:,end));
    sem_NoRibosomes = std(pre_ribosomeDensity(:,end))./ sqrt(nRepetitions);
    
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
    
    %% creating an intensity vector removing values below a detection threshold.
    threshold = 5; % number of
    vec_intensity_Simulation =  intensityVector(:,end);
    vec_intensity_Simulation (vec_intensity_Simulation<threshold)= [];
    vec_intensity_Simulation = vec_intensity_Simulation/10;
    intensity_Simulation{g} = vec_intensity_Simulation;
    
    % normalizing intensity vector
%     for k = 1 : nRepetitions
%         intensityVector(k,:) =  intensityVector(k,:) ./ max(intensityVector(k,:) );
%     end
    
    %% time vector
    time = linspace(0,totalSimulationTime-nonConsideredInitialSimulationTime,totalSimulationTime-nonConsideredInitialSimulationTime);
    maximumNoExperimentalRepetitions=20;
%     [lags, simulation_autocorrelation,simulation_sd_autocorrelation,simulation_sem_autocorrelation] = AutoCorrelation_Simple(intensityVector,maximumNoExperimentalRepetitions);
    [lags, simulation_autocorrelation,simulation_sd_autocorrelation,simulation_sem_autocorrelation] = AutoCorrelation_Simple_normalized_Eq8(intensityVector,maximumNoExperimentalRepetitions);
%     [lags, simulation_autocorrelation,simulation_sd_autocorrelation,simulation_sem_autocorrelation] = AutoCorrelation_Simple_noNormalized(intensityVector,maximumNoExperimentalRepetitions);

    fitAC (g) = computing_Obj_AC (lags,simulation_autocorrelation,simulation_sem_autocorrelation, fileName_AC,maxlag,samplingRateInSeconds, geneFile,folderName, plottingCondition);
    [~,~,fit_KSDistance(g)] = kstest2(intensity_Simulation{g}, renormalized_intensities{g});

    if plottingCondition ==1
        plotCombinationDistribution(geneFile,intensity_Simulation{g},folderName)
    end
end
% fitVal = sum(fit_KSDistance) + sum(fitAC);

fitVal(1) =fit_KSDistance(1);
fitVal(2) =fit_KSDistance(2);
fitVal(3) =fit_KSDistance(3);


fitVal(4) =fitAC(1);
fitVal(5) =fitAC(2);
fitVal(6) =fitAC(3);



% fitVal = round(fitVal,3);
end