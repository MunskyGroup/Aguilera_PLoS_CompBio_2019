function [HT_ke_SR_1s,HT_ke_SR_3s,HT_ke_SR_10s,HT_ke_SR_20s] = function_calculate_HT_ke_REP_2 (parametersModel,tagPositions,vector_NoSpots,meanHT_Delay)
%% Deffing parameter values
nRepetitions = vector_NoSpots(end);
evaluated_no_spots = length(vector_NoSpots);
nonConsideredInitialSimulationTime = 2000;
totalSimulationTime = 6300+ nonConsideredInitialSimulationTime;
timePerturbationApplication= nonConsideredInitialSimulationTime+300;
evaluatingInhibitor = 1; evaluatingFRAP =0;
tP=timePerturbationApplication-nonConsideredInitialSimulationTime;

%% Running Simulations
intensityVector = zeros (nRepetitions,totalSimulationTime-nonConsideredInitialSimulationTime );
parfor k = 1 : nRepetitions
    timePerturbationApplication_withDelay = timePerturbationApplication + randn*meanHT_Delay;
    I =  solve_direct_optimized_2(nonConsideredInitialSimulationTime,totalSimulationTime, parametersModel, timePerturbationApplication_withDelay,evaluatingInhibitor,evaluatingFRAP,tagPositions);
    intensityVector (k,:)= I;
end
geneLength=length(parametersModel)-2;

%% Normalizing intensity vectors
for k = 1 : nRepetitions
    intensityVector(k,:) =  intensityVector(k,:) ./ max(intensityVector(k,:) );
end
%% Prealocating vectors
time_FrameRate_1sec  = linspace(0,totalSimulationTime-nonConsideredInitialSimulationTime,totalSimulationTime-nonConsideredInitialSimulationTime);
reduction_vector=[1,3,10,20];
no_samplingRates=length(reduction_vector);
HT_ke_SR = zeros (1, no_samplingRates);
totalFrames = 300;
totalSimulationTime = reduction_vector*totalFrames + tP;

HT_ke_SR_1s = zeros(1,evaluated_no_spots);
HT_ke_SR_3s = zeros(1,evaluated_no_spots);
HT_ke_SR_10s = zeros(1,evaluated_no_spots);
HT_ke_SR_20s = zeros(1,evaluated_no_spots);


for j=1: evaluated_no_spots
    %% Calculating elongation rates using different frame rates
    for r =1: no_samplingRates
        
        % time reduction
        reduction = reduction_vector(r);
        tp_r= round(tP/reduction);
        time_FrameRate = downsample(time_FrameRate_1sec([1,1:totalSimulationTime(r)-1]),reduction);
        
        % intensity reduction and normalization
        intensityVector_Reduction_noSpots =  intensityVector([1:vector_NoSpots(j)],[1:totalSimulationTime(r)]);
        downsampled_FrameRate_intensityVector = downsample (intensityVector_Reduction_noSpots',reduction)';
        mean_Intensity = mean (downsampled_FrameRate_intensityVector);
        mean_Intensity = mean_Intensity./mean(mean_Intensity([1:tp_r]));
        
        % Linear fit and tau determination
        try
            warning('off','all')
            linreg = 1-((tagPositions(end)/2)/geneLength); %beginning of linear region
            aT1 = find(mean_Intensity(tp_r:end)<linreg,1,'first');
            aT2 = find(mean_Intensity(tp_r:end)<0.05,1,'first');
            
            [p,S] = polyfit(time_FrameRate([tp_r+aT1:1:tp_r+aT2]), mean_Intensity([tp_r+aT1:1:tp_r+aT2]),1);
            [y_fit,~] = polyval(p,time_FrameRate([tp_r+aT1:1:tp_r+aT2+round(100/reduction)]),S);
            
            warning('off','all')
            
            runOff_time= time_FrameRate(y_fit<=0.05)+aT1;
            runOff_time= runOff_time(1);
            HT_ke_SR(1,r) = (geneLength*0.95)./runOff_time;
            
        catch
            % HT_ke_SR(1,r) =nan;
            HT_ke_SR(1,r) =(geneLength*0.95)./time_FrameRate(end);
        end
        
        %         figure(9)
        %          hold on
        %          plot(time_FrameRate,mean_Intensity)
        
        %         figure(10)
        %         plot (time_FrameRate,mean_Intensity,'-k','LineWidth',2)
        %         hold on
        %         plot(time_FrameRate([tp_r+aT1:1:tp_r+aT2+round(100/reduction)]),y_fit,'-r','LineWidth',2)
        %         plot([tp_r,tp_r],[0,1.3],'-b','LineWidth',2)
        %         plot([200,600],[0.05,0.05],'-g','LineWidth',2)
        %         ylim([-0.05,1.2])
        %         xlim([280,550])
        %         xlabel('time','fontsize',20)
        %         ylabel('Intensity','fontsize',20)
        %         title ('HT assays','fontsize',20)
        %
        
    end
    HT_ke_SR_1s(1,j) = HT_ke_SR(1,1);
    HT_ke_SR_3s(1,j) = HT_ke_SR(1,2);
    HT_ke_SR_10s(1,j) = HT_ke_SR(1,3);
    HT_ke_SR_20s(1,j) = HT_ke_SR(1,4);
end

%% making inf values equal to nan
HT_ke_SR_1s(isinf(HT_ke_SR_1s)==1)=0;
HT_ke_SR_3s(isinf(HT_ke_SR_3s)==1)=0;
HT_ke_SR_10s(isinf(HT_ke_SR_10s)==1)=0;
HT_ke_SR_20s(isinf(HT_ke_SR_20s)==1)=0;
end