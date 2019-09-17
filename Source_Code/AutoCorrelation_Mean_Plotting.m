function [lags, mean_autocorrelation] = AutoCorrelation_Mean_Plotting(intensityVector,nRepetitions,samplingRateInSeconds)
nPoints = size(intensityVector,2);
lags = linspace(1,nPoints,nPoints);

% Substracting the mean from each trace
for k = 1 : nRepetitions
    intensityVector(k,:) = (intensityVector(k,:) - mean(intensityVector(k,:)))'/std(intensityVector(k,:));
end

%% Code to Calculate autorocorrelations
for k = 1 : nRepetitions
    autocorrelation = xcorr (intensityVector(k,:),'unbiased');
    autocorrelation = autocorrelation(nPoints:end);
    autocorrelation = autocorrelation/autocorrelation(1);
    preautocorrelation (k,:)  = autocorrelation;
end

%% Normalizing the autocorrelation with respect to the first element.
mean_autocorrelation = mean (preautocorrelation);
mean_autocorrelation = mean_autocorrelation/mean_autocorrelation(1);

% downsampling
lags = downsample(lags,samplingRateInSeconds);
lags(1) =0;
mean_autocorrelation = downsample(mean_autocorrelation, samplingRateInSeconds);

% reporting to the maximum lag
lags (1)=0;

% mean_autocorrelation = mean_autocorrelation(1,1:maxlag);


end