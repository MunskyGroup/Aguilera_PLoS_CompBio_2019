function [FRAP_ke_SR_1s,FRAP_ke_SR_3s,FRAP_ke_SR_10s,FRAP_ke_SR_20s] = function_calculate_FRAP_ke_REP_2 (parametersModel,tagPositions,vector_NoSpots)
nRepetitions = vector_NoSpots(end);
evaluated_no_spots = length(vector_NoSpots);
%% Deffing parameter values
nonConsideredInitialSimulationTime = 2000;
totalSimulationTime = 6300+ nonConsideredInitialSimulationTime;
timePerturbationApplication= nonConsideredInitialSimulationTime+300;
tP = (timePerturbationApplication - nonConsideredInitialSimulationTime);
evaluatingInhibitor = 0;
evaluatingFRAP =1;

%% Creating intensity vectors
intensityVector = zeros (nRepetitions,totalSimulationTime-nonConsideredInitialSimulationTime );
%% Running Simulations
parfor k = 1 : nRepetitions
    I = solve_direct_optimized_FRAP_2(nonConsideredInitialSimulationTime,totalSimulationTime, parametersModel, timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP,tagPositions);
    intensityVector (k,:)= I;
end
%% Normalizing intensity vectors
meanVal =  mean (intensityVector(:,end));
for k = 1 : nRepetitions
    intensityVector (k,:)= intensityVector (k,:)./meanVal;
end
geneLength=length(parametersModel)-2;

%% Prealocating vectors
time_FrameRate_1sec  = linspace(0,totalSimulationTime-nonConsideredInitialSimulationTime,totalSimulationTime-nonConsideredInitialSimulationTime);
reduction_vector=[1,3,10,20];
no_samplingRates=length(reduction_vector);
FRAP_ke_SR = zeros (1, no_samplingRates);

totalFrames = 300;
totalSimulationTime = reduction_vector*totalFrames + tP;

FRAP_ke_SR_1s = zeros(1,evaluated_no_spots);
FRAP_ke_SR_3s = zeros(1,evaluated_no_spots);
FRAP_ke_SR_10s = zeros(1,evaluated_no_spots);
FRAP_ke_SR_20s = zeros(1,evaluated_no_spots);

for j=1: evaluated_no_spots
    %% Calculating intensities for different frame rates
    for r =1: no_samplingRates
        
        % time reduction
        reduction=reduction_vector(r);
        tp_r= round(tP/reduction);
        time_FrameRate = downsample(time_FrameRate_1sec([1,1:totalSimulationTime(r)-1]),reduction);
        
        % intensity reduction and normalization
        intensityVector_Reduction_noSpots =  intensityVector([1:vector_NoSpots(j)],[1:totalSimulationTime(r)]);
        downsampled_FrameRate_intensityVector = downsample (intensityVector_Reduction_noSpots',reduction)';
        mean_Intensity = mean (downsampled_FrameRate_intensityVector);
        
        % Linear fit and tau determination
        try
            warning('off','all')
            
            linreg = ((tagPositions(end)/2)/geneLength); %beginning of lionear region
            
            
            aT1 = find(mean_Intensity(tp_r+1:end)>linreg,1,'first');
            aT2 = find(mean_Intensity(tp_r+1:end)>0.95,1,'first');
            
            [p,S] = polyfit(time_FrameRate([tp_r+aT1:1:tp_r+aT2]), mean_Intensity([tp_r+aT1:1:tp_r+aT2]),1);
            [y_fit,~] = polyval(p,time_FrameRate([tp_r+aT1:1:end]),S);
            warning('off','all')
            
            recovery_time= time_FrameRate(y_fit>0.95)+aT1;
            recovery_time = recovery_time(1);
            FRAP_ke_SR(1,r)  = (geneLength*0.95)./recovery_time;
            
        catch
            FRAP_ke_SR(1,r) =(geneLength*0.95)./time_FrameRate(end);
        end
        
        %                 figure(9)
        %                 hold on
        %                 plot(time_FrameRate,mean_Intensity)
        %
        %                         figure(10)
        %                         plot (time_FrameRate,mean_Intensity,'-k','LineWidth',2)
        %                         hold on
        %                         plot(time_FrameRate([tp_r+aT1:1:end]),y_fit,'-r','LineWidth',2)
        %                         plot([tp_r,tp_r],[0,1.3],'-b','LineWidth',2)
        %                         plot([200,600],[0.95,0.95],'-g','LineWidth',2)
        %                         ylim([-0.1,1.2])
        %                         xlim([280,550])
        %                         xlabel('time','fontsize',20)
        %                         ylabel('Intensity','fontsize',20)
        %                         title ('FRAP assays','fontsize',20)
    end
    FRAP_ke_SR_1s(1,j) = FRAP_ke_SR(1,1);
    FRAP_ke_SR_3s(1,j) = FRAP_ke_SR(1,2);
    FRAP_ke_SR_10s(1,j) = FRAP_ke_SR(1,3);
    FRAP_ke_SR_20s(1,j) = FRAP_ke_SR(1,4);
end

%% making nan inf numbers
FRAP_ke_SR_1s(isinf(FRAP_ke_SR_1s)==1)=0;
FRAP_ke_SR_3s(isinf(FRAP_ke_SR_3s)==1)=0;
FRAP_ke_SR_10s(isinf(FRAP_ke_SR_10s)==1)=0;
FRAP_ke_SR_20s(isinf(FRAP_ke_SR_20s)==1)=0;
end
