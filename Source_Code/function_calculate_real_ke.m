function [real_ke,associaiton_Time] = function_calculate_real_ke (parametersModel)
% Deffing parameter values
nonConsideredInitialSimulationTime = 1000;
totalSimulationTime = 5000+ nonConsideredInitialSimulationTime;
 geneLength =length(parametersModel)-2;


timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;
%% Running Simulations
[~,~,RibosomePositions] = solve_direct_optimized(totalSimulationTime, parametersModel, timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP);
%real_ke = realElongationRates(RibosomePositions,nonConsideredInitialSimulationTime,geneLength);
[real_ke,associaiton_Time] = realElongationRates(RibosomePositions,geneLength);

end