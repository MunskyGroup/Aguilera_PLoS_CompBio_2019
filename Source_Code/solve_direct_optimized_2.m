function [intensityVector] = solve_direct_optimized_2(nonConsideredInitialSimulationTime,nPoints, k,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP,tagPositions)

TimeVectorFixedSize = linspace (0,nPoints,nPoints);
geneLength = length(k) - 2;
%% Runnin the SSA
try
RibosomePositions = SSA_longNames_mex(k,TimeVectorFixedSize,timePerturbationApplication,evaluatingInhibitor,0);
catch
RibosomePositions = SSA_longNames(k,TimeVectorFixedSize,timePerturbationApplication,evaluatingInhibitor,0);
end
% Remove rows with only zeros
RibosomePositions = RibosomePositions(any(RibosomePositions,2),:);
%%% In case of using FRAP. Remove ribosomes that initated before the timePerturbationApplication
if evaluatingFRAP ==1
    InitialRib = nnz( RibosomePositions(:,timePerturbationApplication));
    nRib =InitialRib;
    for tp=timePerturbationApplication:nPoints-1
        if nRib>0
            if any (RibosomePositions(1:nRib,tp) > RibosomePositions(1:nRib,tp+1))
                nRib = nRib-1;
            end
            if nRib >0
                RibosomePositions(1:nRib,tp) =0;
            end
        end
    end
end
%% Reconstructing the SSA output: from ribosome dynamics to codon occupancy
sout_temp = zeros(length(TimeVectorFixedSize) , geneLength);
for i = 1:length(TimeVectorFixedSize)
    v = RibosomePositions(RibosomePositions(:,i) ~= 0,i);
    if ~isempty(v)
        sout_temp(i,v) = 1;
    end
end
% s_out = sout_temp;
% t_out = TimeVectorFixedSize;

s_out_Short = sout_temp(nonConsideredInitialSimulationTime:nPoints,:);

probePositionVector = zeros (1, geneLength);
for i =1: length(tagPositions)
    probePositionVector(1, tagPositions(i):end) = 1*i;
end
I = zeros (1,nPoints-nonConsideredInitialSimulationTime);
for tp =1: nPoints-nonConsideredInitialSimulationTime
    temp_simulationOutPut = s_out_Short(tp,1:geneLength);
    I(tp)= sum (probePositionVector .* temp_simulationOutPut );
end
intensityVector = I;
end
