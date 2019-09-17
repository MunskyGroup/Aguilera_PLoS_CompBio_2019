function [lags, mean_autocorrelation,SEM] = AutoCorrelation_Simple_2(intensityVector,nRepetitions)
nPoints = size(intensityVector,2);
lags = linspace(1,nPoints,nPoints);
% lags = downsample(lags,samplingRateInSeconds);

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
G0 = mean_autocorrelation(1);
mean_autocorrelation = mean_autocorrelation/G0;
% Calculating SEM and SD
SD = std (preautocorrelation)/ G0^2;
SEM = SD/sqrt(nRepetitions);
% norm_mean_autocorrelation(isnan(norm_mean_autocorrelation))=0;
SEM(isnan(SEM))=0;

end