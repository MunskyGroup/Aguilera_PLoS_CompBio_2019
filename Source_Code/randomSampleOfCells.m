function [randIntensities,AllIntensity] = randomSampleOfCells(g,nRep)
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

for j=1:nRep
    NN = ceil(length(StrData)*rand(1,floor(length(StrData))));
    StrData_New = StrData(NN);
    %% calculating mean intensities in arbitraty units
    rand_vec_Intensity = [];
    for k =1:length(StrData_New)
        rand_vec_Intensity = [rand_vec_Intensity,StrData_New(1,k).trajectory'];
    end
    randIntensities{j} = rand_vec_Intensity;
end

%% calculating all intensities
AllIntensity =[];
i=1;
while size(find(rawData(:,1)==i))>0
    AllIntensity =  [AllIntensity,StrData(1,i).trajectory'];
        i=i+1;
end

end