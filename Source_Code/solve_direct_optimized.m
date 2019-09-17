function [t_out,s_out,RibosomePositions,s_out_neg] = solve_direct_optimized(nPoints, k,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP)
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
RibosomePositions_neg = 0*RibosomePositions;
if evaluatingFRAP ==1
    InitialRib = nnz( RibosomePositions(:,timePerturbationApplication));
    nRib =InitialRib;
    Ribs_at_FRAP = RibosomePositions(1:nRib,timePerturbationApplication);
    for tp=timePerturbationApplication:nPoints-1
        if nRib>0
            if any (RibosomePositions(1:nRib,tp) > RibosomePositions(1:nRib,tp+1))
                nRib = nRib-1;
            end
            if nRib > 0
                %                 RibosomePositions(1:nRib,tp) =0;
                RibosomePositions_neg(1:nRib,tp) = Ribs_at_FRAP(length(Ribs_at_FRAP)-nRib+1:length(Ribs_at_FRAP),1);
            end
        end
    end
end
%% Reconstructing the SSA output: from ribosome dynamics to codon occupancy
sout_temp = zeros(length(TimeVectorFixedSize) , geneLength);
sout_temp_neg = zeros(length(TimeVectorFixedSize) , geneLength);
for i = 1:length(TimeVectorFixedSize)
    v = RibosomePositions(RibosomePositions(:,i) ~= 0,i);
    v_neg = RibosomePositions_neg(RibosomePositions_neg(:,i) ~= 0,i);
    if ~isempty(v)
        sout_temp(i,v) = 1;
    end
    if ~isempty(v_neg)
        sout_temp_neg(i,v_neg) = 1;
    end
end
s_out = sout_temp;
s_out_neg = sout_temp_neg;
t_out = TimeVectorFixedSize;
