function [SEM] = AutoCorrelation_SEM_Plotting(intensityVector,nExpRepetitions,samplingRateInSeconds)


downsampled_intensityVector = downsample(intensityVector',samplingRateInSeconds);
downsampled_intensityVector=downsampled_intensityVector';
nPoints = size(downsampled_intensityVector,2);


lags = linspace(1,nPoints,nPoints);
% Substracting the mean from each trace
for k = 1 : nExpRepetitions
    downsampled_intensityVector(k,:) = (downsampled_intensityVector(k,:) - mean(downsampled_intensityVector(k,:)))'/std(intensityVector(k,:));
end

%% Code to Calculate autorocorrelations
for k = 1 : nExpRepetitions
    autocorrelation = xcorr (downsampled_intensityVector(k,:),'unbiased');
    autocorrelation = autocorrelation(nPoints:end);
    autocorrelation = autocorrelation/autocorrelation(1);
    matrixAutocorrelation (k,:)  = autocorrelation;
end

%% Normalizing the autocorrelation with respect to the first element.
% Calculating the mean form all repetitions
mean_autocorrelation = mean (matrixAutocorrelation);
mean_autocorrelation = mean_autocorrelation/mean_autocorrelation(1);

% Calculating SEM and SD
SEM = std(matrixAutocorrelation)./sqrt(nExpRepetitions);
SEM= SEM/(mean_autocorrelation(1)^2);

% size(mean_autocorrelation)
% size(SEM)


% norm_mean_autocorrelation(isnan(norm_mean_autocorrelation))=0;
SEM(isnan(SEM))=0;

end