function [lags_exp, mean_autocorrelation, covariance_Autocorrelation, sd_autocorrelation,nTraces,matrix_Autocorrelation] = autoCorrelation_Data_matrices(fileName, samplingRateInSeconds, maxlag)

[~,sheet_name]=xlsfinfo(fileName);
for i=1:numel(sheet_name)
    %% Extracting data and removing zeros
    data{i} = xlsread(fileName,sheet_name{i});
    temp1 = data{i}(:,2);
    StrData(1,i).trajectory = temp1;
end

% IKeep = zeros(1,nTraces,'logical');
% for i=1:nTraces
%     if length(StrData(i).trajectory)>=50
%         IKeep(i)=1;
%     end
% end
% StrData = StrData(IKeep);
nTraces = size(StrData,2);


%% Calculating AC ignoring zeros
removingZeros =0; % 0 keeps zeros ; 1 removes zeros

% The experimental data may contain unreadable intensities. In our case these are represented by zeros.
% to solve this problem one option is to remove zeros form the dataset and lump adjacent time points.
% the second option is to keep this zeros as part of the calculation.
% Unfortunatly both approaches induce some error and affect the autocovariance calculation.
% We tested both approaches and similar results are obtained.

if removingZeros ==0  % AC Calculation keeping zeros.
    for i=1:nTraces
        meanInt(i)= mean(StrData(i).trajectory);
        stdInt(i)= std(StrData(i).trajectory);
    end
    
    for j =1: nTraces
        nPoints = length (StrData(j).trajectory);
        autocorrelation =[];
        % loading experimental data and substracting the mean.
        IntensityNormalized = (StrData(j).trajectory - meanInt(j))'/stdInt(j);
        autocorrelation = xcorr (IntensityNormalized,'unbiased');
        autocorrelation = autocorrelation(nPoints:end);
        if length (autocorrelation)>= maxlag
            autocorrelation = autocorrelation(1:maxlag);
        else
            autocorrelation(end:maxlag) = 0; % filling last elements with zeros
        end
        matrix_Autocorrelation(j,:)= autocorrelation;
    end
    mean_autocorrelation= mean(matrix_Autocorrelation);
    
else     % AC Calculation removing zeros and lumping values.
    
    for i=1:nTraces
        temp =[];
        temp =StrData(i).trajectory;
        temp(temp==0)=NaN;
        meanInt(i)= nanmean(temp);
        stdInt(i)= nanstd(temp);
    end
    
    for j =1: nTraces
        autocorrelation =[];
        tempInt =[];
        tempInt= StrData(j).trajectory ;
        tempInt(tempInt==0)=[];
        nPoints = length (tempInt);
        % loading experimental data and substracting the mean.
        IntensityNormalized = (tempInt- meanInt(j))'/stdInt(j);
        autocorrelation = xcorr (IntensityNormalized,'unbiased');
        autocorrelation = autocorrelation(nPoints:end);
        if length (autocorrelation)>= maxlag
            autocorrelation = autocorrelation(1:maxlag);
        else
            autocorrelation(end:maxlag) = 0; % filling last elements with zeros
        end
        matrix_Autocorrelation(j,:)= autocorrelation;
    end
    mean_autocorrelation= mean(matrix_Autocorrelation);
end


X = [2 3 4 5]';
if length(X)==1
    G0 = mean_autocorrelation(X+1);
else
    V = mean_autocorrelation(X+1)';
    G0 = interp1(X,V,0,'linear','extrap');
end
mean_autocorrelation = mean_autocorrelation/G0;
mean_autocorrelation(1) = 1;

sd_autocorrelation = std(matrix_Autocorrelation)/sqrt(nTraces);
sd_autocorrelation = sd_autocorrelation/ G0^2;
sd_autocorrelation(1) = 0;

lags_exp = linspace(1,maxlag,maxlag)*samplingRateInSeconds;
lags_exp(1) = 0;

% matrix_Autocorrelation(:,1)=[]; % removing the first element
covariance_Autocorrelation = cov(matrix_Autocorrelation);
covariance_Autocorrelation = covariance_Autocorrelation/G0^2;
covariance_Autocorrelation(1,1)=0;

end