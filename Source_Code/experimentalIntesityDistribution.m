%function [renormalized_intensities] = experimentalIntesityDistribution 
function [experimentalData] = experimentalIntesityDistribution 
% returns a vector with the intensity for the experimental data
for g = 1:3

    if g==1
filename = 'intensities_KDM5B.xlsx'; rawData = xlsread(filename,1,'A1:B6435'); scalingFactor = 4.99;% KDM5B
    end
    if g==2
 filename = 'intensities_Bact.xlsx'; rawData = xlsread(filename,1,'A1:B3973'); scalingFactor=3; % B-actin
    end
    if g==3
filename = 'intensities_H2B.xlsx'; rawData = xlsread(filename,1,'A1:B751'); scalingFactor=2.12; % H2B
    end

%% extracting data from excell files
i=1;
while size(find(rawData(:,1)==i))>0
    %% Creating structure of data.
    StrData(1,i).trajectory = (rawData(rawData(1:end,1) == i, 2));
       i=i+1;
end

%% calculating mean intensities in arbitraty units
AllIntensity =[];
N_frames_for_cell_dists = 50;
for k =1:N_frames_for_cell_dists
% rint = ceil(rand*length(StrData(1,k).trajectory));
% AllIntensity(k) =  StrData(1,k).trajectory(rint);
AllIntensity =  [AllIntensity,StrData(1,k).trajectory'];
end
experimentalData{g} = AllIntensity;
end
end