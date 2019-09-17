clear all; close all; clc;
cd Source_Code

%% Code to calculate parameter values between the different model implementations.
%% This code reproduces the values given in table 1.
%% All results are directly printed to screen.

%% Deffing parameter values
nRepetitions = 500;
nonConsideredInitialSimulationTime = 4000;
totalSimulationTime = 4000+ nonConsideredInitialSimulationTime;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;

for rep=1:3
    
    for g =1:3
        
        %% loading experimental data
        if g==1
            fileName = 'KDM5B.xls'; samplingRateInSeconds= 10; samplingForPlotting=1; maxlag =1000;geneFile = 'KDM5B_withTags.txt';bar_ke = 10; ki=0.03;   fileName2 = 'KDM5B_AuCo.xls';
        elseif g==2
            fileName = 'bact.xls'; samplingRateInSeconds= 3; samplingForPlotting=1; maxlag =300;  geneFile = 'Bactin_withTags.txt'; bar_ke = 10;  ki=0.03; fileName2 = 'ActB_AuCo.xls';
        elseif g==3
            fileName = 'H2B.xls'; samplingRateInSeconds= 1; samplingForPlotting=1; maxlag =100;  geneFile = 'H2B_withTags.txt'; bar_ke = 10; ki=0.03; fileName2 = 'H2B_AuCo.xls';
        end
        
        %% Creating Model from Gene sequence
        [generated_geneSequence, typeOfTag,tagPositions,parametersModel ,geneLength] = sequenceAnalyzer(geneFile,ki,bar_ke);
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
        
        intensity_UMP =[];
        intensity_UMP=  intensityVector/10;
        mean_Intensity_UMP(g) = mean (mean (intensity_UMP,2));
        var_Intensity_UMP(g) = mean ((std (intensity_UMP,0,2)).^2);

        
        [real_ke,associaiton_Time] = function_calculate_real_ke (parametersModel);
        tau_period (g) = associaiton_Time;
    end
    
    %%
    meanInt (rep,:)= mean_Intensity_UMP;
    var_Intensit (rep,:) = var_Intensity_UMP;
    tau (rep,:) = tau_period;
    
end

mean_Int_final = mean(meanInt)
std_Int_final = std(meanInt)

mean_var_Intensit = mean(var_Intensit)
std_var_Intensit = std(var_Intensit)

mean_tau_final = mean(tau)
std_tau_final = std(tau)


%% loading experimental data
bar_ke = 10; ki=0.03;  

fileName = 'KDM5B.xls'; samplingRateInSeconds= 10; samplingForPlotting=1; maxlag =1000;geneFile = 'KDM5B_withTags.txt'; bar_ke = 10;
[generated_geneSequence, typeOfTag,tagPositions,parametersModel ,geneLength_kdm5b] = sequenceAnalyzer(geneFile,ki,bar_ke);
[codonUsage ] = sequenceUsage(geneFile);
tau_kdm5b = 1/bar_ke * codonUsage;

fileName = 'bact.xls'; samplingRateInSeconds= 3; samplingForPlotting=1; maxlag =300;  geneFile = 'Bactin_withTags.txt'; bar_ke = 10;
[generated_geneSequence, typeOfTag,tagPositions,parametersModel ,geneLength_bact] = sequenceAnalyzer(geneFile,ki,bar_ke);
[codonUsage ] = sequenceUsage(geneFile);
tau_bac = 1/bar_ke * codonUsage;

fileName = 'H2B.xls'; samplingRateInSeconds= 1; samplingForPlotting=1; maxlag =100;  geneFile = 'H2B_withTags.txt'; bar_ke = 10;
[generated_geneSequence, typeOfTag,tagPositions,parametersModel ,geneLength_h2b] = sequenceAnalyzer(geneFile,ki,bar_ke);
[codonUsage ] = sequenceUsage(geneFile);
tau_h2b = 1/bar_ke * codonUsage;

% % calculating mu_Int corrected by mean probe position

% mu_Int_kdm5b = tau_kdm5b * ki* ( 1 - (mean (tagPositions) / geneLength_kdm5b))
% var_Int_kdm5b = tau_kdm5b * ki* ( 1 - (mean (tagPositions) / geneLength_kdm5b)).^2
% 
% 
% mu_Int_bact = tau_bac * ki*  ( 1 - (mean (tagPositions) / geneLength_bact))
% var_Int_bact = tau_bac * ki*  ( 1 - (mean (tagPositions) / geneLength_bact)).^2
% 
% 
% mu_Int_h2b = tau_h2b * ki*  ( 1 - (mean (tagPositions) / geneLength_h2b))
% var_Int_h2b = tau_h2b * ki*  ( 1 - (mean (tagPositions) / geneLength_h2b)).^2

%% Second for of correction
clc

mu_Int_kdm5b = tau_kdm5b * ki* ( 1 - (tagPositions(end) / (2*geneLength_kdm5b)))
var_Int_kdm5b = tau_kdm5b * ki* ( 1 - ((2/3) *tagPositions(end)/ geneLength_kdm5b))

mu_Int_bact = tau_bac * ki*  ( 1 - (tagPositions(end)/ (2*geneLength_bact)))
var_Int_bact = tau_bac * ki*  ( 1 - ((2/3) *tagPositions(end) / geneLength_bact))

mu_Int_h2b = tau_h2b * ki*  ( 1 - ( tagPositions(end)/ (2*geneLength_h2b)))
var_Int_h2b = tau_h2b * ki*  ( 1 - ((2/3) * tagPositions(end)/ geneLength_h2b))

cd ..

