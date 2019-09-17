
function [real_ke,associaiton_Time] = realElongationRates(RibosomePositions,geneLength)
% Converting Zeros to Nan in Ribosome Positions
geneLength = geneLength;
%nonConsideredInitialSimulationTime =1000;

RibosomePositionsNaN= RibosomePositions;
RibosomePositionsNaN(RibosomePositionsNaN==0) =nan;
%RibosomePositionsNaN= RibosomePositionsNaN(:,[nonConsideredInitialSimulationTime+1:end]);

tp = size(RibosomePositionsNaN,2);
intiation = zeros(1,tp);
termination = zeros(1,tp);

for tp = 1 : tp-1
    if isnan(min(RibosomePositionsNaN(:,tp)))==1 && isnan(min(RibosomePositionsNaN(:,tp+1)))==0
        intiation(1,tp) = 1;
    else
        intiation(1,tp) = min(RibosomePositionsNaN(:,tp)) > min(RibosomePositionsNaN(:,tp+1));
    end
    
    if isnan(max(RibosomePositionsNaN(:,tp)))==0 && isnan(max(RibosomePositionsNaN(:,tp+1)))==1
        termination(1,tp) =1;
    else
        termination(1,tp) = max(RibosomePositionsNaN(:,tp)) > max(RibosomePositionsNaN(:,tp+1));
    end
end

%% Calculating Real Elongation Rate
tp = size(RibosomePositionsNaN,2);
counter_1=1;
elongationTime =0;
for tp = 1 : tp-1
    if intiation(1,tp) ==1
        trackedRibosome= min(RibosomePositionsNaN(:,tp+1));
        tempValue_Ribosome = trackedRibosome;
        counter_2 =2;
        terminationCondition=0;
        try
        while trackedRibosome < geneLength && terminationCondition==0
            Index = find(RibosomePositionsNaN(:,tp+counter_2)>= trackedRibosome);
            if isempty(Index)==0
            tempValue_Ribosome = RibosomePositionsNaN(Index,tp+counter_2);
            tempValue_Ribosome = sort(tempValue_Ribosome);
            trackedRibosome = tempValue_Ribosome(1);
            counter_2 = counter_2 +1;
            else
                terminationCondition =1;
            end
        end
        elongationTime(1,counter_1) =counter_2;
        counter_1=counter_1+1;
        catch
        end
    end
end


if elongationTime ~=0
real_ke = round(geneLength/mean(elongationTime),1);
associaiton_Time = mean(elongationTime);
else
    real_ke = 0;
associaiton_Time = inf;
end
end