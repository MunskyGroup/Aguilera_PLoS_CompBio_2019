function [lags, mean_autocorrelation,SEM] = AutoCorrelation_Simple(intensityVector,nRepetitions,samplingRateInSeconds)
nPoints = size(intensityVector,2);
lags = linspace(1,nPoints,nPoints)*samplingRateInSeconds;
lags = downsample(lags,samplingRateInSeconds);

lags (1)= 0;

% Substracting the mean from each trace
for k = 1 : nRepetitions
    intensityVector(k,:) = (intensityVector(k,:) - mean(intensityVector(k,:)))'/std(intensityVector(k,:));
end

%% Code to Calculate autorocorrelations
for k = 1 : nRepetitions
    autocorrelation = xcorr (intensityVector(k,:),'unbiased');
    autocorrelation = autocorrelation(nPoints:end);
%     autocorrelation = autocorrelation/autocorrelation(1);    
    preautocorrelation (k,:)  = autocorrelation;
end

%% Normalizing the autocorrelation with respect to the first element.
% Calculating the mean form all repetitions
mean_autocorrelation = mean (preautocorrelation);
mean_autocorrelation = mean_autocorrelation/mean_autocorrelation(1);
mean_autocorrelation = downsample(mean_autocorrelation, samplingRateInSeconds);

% Calculating SEM and SD
SEM = (std(preautocorrelation)./sqrt(nRepetitions))/mean_autocorrelation(1);
SEM = downsample(SEM, samplingRateInSeconds);

% norm_mean_autocorrelation(isnan(norm_mean_autocorrelation))=0;
SEM(isnan(SEM))=0;

end