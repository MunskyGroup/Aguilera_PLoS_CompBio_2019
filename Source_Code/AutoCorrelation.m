function [lags, norm_mean_autocorrelation,normalized_sem] = AutoCorrelation(intensityVector,nRepetitions,tAC,reduction)
nPoints = size(intensityVector,2);
lags = linspace(1,nPoints,nPoints);
for k = 1 : nRepetitions
    intensityVector(k,:) = (intensityVector(k,:) - mean(intensityVector(k,:)))'/std(intensityVector(k,:)) ;
end
for k = 1 : nRepetitions
    autocorrelation = xcorr (intensityVector(k,:),'unbiased');
    autocorrelation = autocorrelation(nPoints:end);
    preautocorrelation (k,:)  = autocorrelation;
end
mean_autocorrelation = mean (preautocorrelation);
mean_autocorrelation = mean_autocorrelation + abs(mean (mean_autocorrelation([tAC:tAC+round(100/reduction)])));
norm_mean_autocorrelation = mean_autocorrelation./mean_autocorrelation(1);
normalized_sem = std(preautocorrelation)/mean_autocorrelation(1);
%% removing Nan
norm_mean_autocorrelation(isnan(norm_mean_autocorrelation))=0;
normalized_sem(isnan(normalized_sem))=0;
end