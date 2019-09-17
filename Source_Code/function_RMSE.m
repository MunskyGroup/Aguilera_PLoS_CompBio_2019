function [RMSE_long,RMSE_medium,RMSE_short] = function_RMSE (calculated_ke,N_Rep_Spots,numElementsLibrary,geneLength,asuumed_ke)
% NoRepetitions = size(calculated_ke,2);
% calculated_ke (isnan(calculated_ke)==1)=0;
for i=1:N_Rep_Spots
    if isnan(calculated_ke(i))==1
calculated_ke(i)= geneLength(i)/300;
    end
end

%% prealocating memory
 RMSE_short = zeros(1,N_Rep_Spots);
 RMSE_medium = zeros(1,N_Rep_Spots);
 RMSE_long = zeros(1,N_Rep_Spots);

%% Calculating correlation coefficients between real elongation and calculated method
for Rep =1: N_Rep_Spots
    %% Deffining counters
    c1=1; c2=1; c3=1;    
    % Separating data sets accorging to length.
    for i =1:numElementsLibrary
        % Select medium size genes
        if geneLength(1,i)>500 && geneLength(1,i)<= 1000 %&& calculated_ke(i,Rep) >0
            ke_medium(1,c2) = calculated_ke(i,Rep);
            c2=c2+1;
            % Select long genes
        elseif geneLength(1,i) >1000 %&& calculated_ke(i,Rep) >0
            ke_long(1,c3) = calculated_ke(i,Rep);
            c3=c3+1;
            % Select small genes
        elseif geneLength(1,i)<=500 %&& calculated_ke(i,Rep) >0
            ke_short(1,c1) = calculated_ke(i,Rep);
            c1=c1+1;
        end
    end
    %% Calculating the Root mean square error
    RMSE_short  (1,Rep) = sqrt(nanmean((asuumed_ke - ke_short).^2));  % 
    RMSE_medium (1,Rep) = sqrt(nanmean((asuumed_ke - ke_medium).^2));  % 
    RMSE_long (1,Rep)   = sqrt(nanmean((asuumed_ke - ke_long).^2));  %
end
% end