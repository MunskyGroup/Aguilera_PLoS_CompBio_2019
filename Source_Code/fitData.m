%% This code is inteded to simulate single molecule stochastic dynamics.
function fitVal = fitData(x,experimentalData, plottingCondition,folderName)

% To facilitate the optimization process the parameters are scaled to
% integer values. Elongation rates are scaled as ke_scale = ke_real *10
% Initiation rates are scaled as ki_scaled = ki_real*1e4;
% This section convert the scaled parameters to real values
k_elongation1 = x(1)/10;
k_initiation_1 = x(2)/1000;
k_initiation_2 = x(3)/1000;
k_initiation_3 = x(4)/1000;

%% Deffing parameter values
nRepetitions = 100;
nonConsideredInitialSimulationTime = 1000;
totalSimulationTime = 3000+ nonConsideredInitialSimulationTime;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;
%% loading experimental data
for g=1:3
    %     c=c-1;
    if g==1
        fileName = 'kdmIntensityData.xls'; samplingRateInSeconds= 10; maxlag =250;geneFile = 'KDM5B_withTags.txt';  ExpRib=4.99; k_initiation =k_initiation_1; k_elongationMean=k_elongation1;nExpRepetitions = 28;
        
    end
    if g==2
        fileName = 'actBIntensityData.xls'; samplingRateInSeconds= 3; maxlag =250;  geneFile = 'Bactin_withTags.txt';   ExpRib=3.01; k_initiation =k_initiation_2;k_elongationMean=k_elongation1;nExpRepetitions = 17;
    end
    if g==3
        fileName = 'h2bIntensityData.xls'; samplingRateInSeconds= 1; maxlag =250;  geneFile = 'H2B_withTags.txt';  ExpRib=2.12; k_initiation=k_initiation_3;k_elongationMean=k_elongation1;nExpRepetitions = 9;
    end
    
    %% Creating Model from Gene sequence
    [~, ~,tagPositions,parametersModel ,geneLength] = sequenceAnalyzer(geneFile,k_initiation,k_elongationMean);
    geneLength=geneLength-1;
    
    %% Stochastic simulations
    parfor k = 1 : nRepetitions
        [~,s_out,~,~] = solve_direct_optimized(totalSimulationTime, parametersModel,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP);
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
    
    %% creating an intensity vector removing values below a detection threshold.
    threshold = 5; % number of
    vec_intensity_Simulation =  intensityVector(:,end);
    vec_intensity_Simulation (vec_intensity_Simulation<threshold)= [];
    vec_intensity_Simulation = vec_intensity_Simulation/10;
    intensity_Simulation{g} = vec_intensity_Simulation;
    
    %% Autocorrelation
    time = linspace(0,totalSimulationTime-nonConsideredInitialSimulationTime,totalSimulationTime-nonConsideredInitialSimulationTime);
    if plottingCondition ==0
        [sim_lags, simulation_autocorrelation,simulation_sem_autocorrelation] = AutoCorrelation_Simple(intensityVector,nRepetitions,samplingRateInSeconds);
    else
        [sim_lags,simulation_autocorrelation] = AutoCorrelation_Mean_Plotting(intensityVector,nExpRepetitions,samplingRateInSeconds);
        [simulation_sem_autocorrelation] = AutoCorrelation_SEM_Plotting(intensityVector,nExpRepetitions,samplingRateInSeconds);
    end
    fit_AC (g) = computing_Obj_AC (sim_lags,simulation_autocorrelation,simulation_sem_autocorrelation, fileName,maxlag,samplingRateInSeconds, geneFile,folderName, plottingCondition);
    
    %% Distributions
    edges =  [0:1:20];
    [hist_sim,~] = histcounts(intensity_Simulation{g},edges, 'Normalization', 'probability');
    [hist_exp,~] = histcounts(experimentalData{g},edges); %;, 'Normalization', 'probability');
    
    % normalize histogram to size of 100 spots per gene (to account for
    % different numbers of spots seen for each gene. This equalizes the
    % emphasis on KDM5B compared to H2B.
    N_spots_used_for_distributions = 100;
    hist_exp = hist_exp*N_spots_used_for_distributions/sum(hist_exp);
    
    hist_sim(hist_sim==0)=1e-10;
    fit_Distributions(g) =  (-  dot(hist_exp,log(hist_sim))) ;
    
    if plottingCondition ==1
        plotCombinationDistribution(geneFile,intensity_Simulation{g},folderName)
        plotTimeCourses (time, intensityVector,nExpRepetitions,geneFile,folderName,samplingRateInSeconds)
        plotTimeCoursesExp (fileName,maxlag,nExpRepetitions,geneFile,folderName,samplingRateInSeconds)
    end
end

% fit_Distributions
% fit_AC
fitVal= sum(fit_Distributions) + sum(fit_AC);
end