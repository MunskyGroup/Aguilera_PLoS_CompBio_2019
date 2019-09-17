function [ac_ke_SR1,ac_ke_SR3,ac_ke_SR10,ac_ke_SR20] = function_calculate_ac_ke_REP_2 (parametersModel,tagPositions,vector_NoSpots)
%% Deffing parameter values
nRepetitions = vector_NoSpots(end);
evaluated_no_spots = length(vector_NoSpots);
nonConsideredInitialSimulationTime = 2000;
totalSimulationTimewithBurning = 6000+ nonConsideredInitialSimulationTime;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;

%% Running Simulations
intensityVector = zeros (nRepetitions,totalSimulationTimewithBurning-nonConsideredInitialSimulationTime );
parfor k = 1 : nRepetitions
    I = solve_direct_optimized_2(nonConsideredInitialSimulationTime,totalSimulationTimewithBurning, parametersModel, timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP,tagPositions);
    intensityVector (k,:)= I;
end
geneLength=length(parametersModel)-2;

%% Prealocating vectors
time_FrameRate_1sec = linspace(0,totalSimulationTimewithBurning-nonConsideredInitialSimulationTime,totalSimulationTimewithBurning-nonConsideredInitialSimulationTime);
reduction_vector=[1,3,10,20];
no_samplingRates=length(reduction_vector);
ac_ke_SR = zeros (1, no_samplingRates);
totalFrames =300;
totalSimulationTime = reduction_vector*totalFrames;

    ac_ke_SR1 = zeros(1,evaluated_no_spots);
    ac_ke_SR3 = zeros(1,evaluated_no_spots);
    ac_ke_SR10 = zeros(1,evaluated_no_spots);
    ac_ke_SR20 = zeros(1,evaluated_no_spots);

%% Calculating elongation rates using different frame rates
for j=1: evaluated_no_spots
    for r =1: no_samplingRates
        % time reduction
        reduction=reduction_vector(r);
        time_FrameRate = downsample(time_FrameRate_1sec([1,1:totalSimulationTime(r)]),reduction);
        % intensity reduction
        intensityVector_Reduction_noSpots =  intensityVector([1:vector_NoSpots(j)],[1:totalSimulationTime(r)]);
        intensityVector_SR = downsample(intensityVector_Reduction_noSpots',reduction)';
        try
            warning('off','all')
            k_elongationMean= mean(parametersModel(2:end-1));
            tAC = round(((geneLength/ k_elongationMean) )/reduction);
            [~, simulation_autocorrelation,~] = AutoCorrelation(intensityVector_SR,vector_NoSpots(j),tAC,reduction);
            Index_SR=find(simulation_autocorrelation <0.05, 1, 'first');
            dwell_time = time_FrameRate (Index_SR);
            warning('off','all')            
            ac_ke_SR(1,r) = (geneLength*0.95)./dwell_time;
            
        catch
            ac_ke_SR(1,r) =(geneLength*0.95)./time_FrameRate(end);
        end
    end
    ac_ke_SR1(1,j) = ac_ke_SR(1,1);
    ac_ke_SR3(1,j) = ac_ke_SR(1,2);
    ac_ke_SR10(1,j) = ac_ke_SR(1,3);
    ac_ke_SR20(1,j) = ac_ke_SR(1,4);
end

ac_ke_SR1(isinf(ac_ke_SR1)==1)=0;
ac_ke_SR3(isinf(ac_ke_SR3)==1)=0;
ac_ke_SR10(isinf(ac_ke_SR10)==1)=0;
ac_ke_SR20(isinf(ac_ke_SR20)==1)=0;
end