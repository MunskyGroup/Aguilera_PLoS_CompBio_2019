function [t_out,s_out,collisions,X_array_RibosomeDynamcis] = solve_direct_optimized_collisions(nPoints, k,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP)


TimeVectorFixedSize = linspace (0,nPoints,nPoints);
 geneLength = length(k) - 2;
 %% Runnin the SSA
    [X_array_RibosomeDynamcis,collisions] = SSA_longNames_collisions(k,TimeVectorFixedSize,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP);
 
 %% Reconstructing the SSA output: from ribosome dynamics to codon occupancy       
    sout_temp = zeros(length(TimeVectorFixedSize) , geneLength);   
    for i = 1:length(TimeVectorFixedSize)
        v = X_array_RibosomeDynamcis(X_array_RibosomeDynamcis(:,i) ~= 0,i);
        if ~isempty(v)
            sout_temp(i,v) = 1;
        end
    end
s_out = sout_temp;
t_out = TimeVectorFixedSize;
