function fitAC = computing_Obj_AC (lags, simulation_autocorrelation,simulation_sem_autocorrelation,fileName,maxlag,samplingRateInSeconds,geneFile,folderName,plottingCondition)

%% Loading experimental data
[experimental_lags,experimental_autocorrelation,covariance_Autocorrelation,experimental_sem,nTraces,matrix_Autocorrelation] = autoCorrelation_Data_matrices(fileName, samplingRateInSeconds, maxlag);

%% Plotting 1 - Autocorrelation.
if plottingCondition==1
    figure('visible', 'off');
    fig1= gcf;
    fig1.PaperUnits = 'inches';
    fig1.PaperPosition = [0, 0, 2, 1.1];
    hold on
    h=errorbar(experimental_lags(1:maxlag), experimental_autocorrelation(1:maxlag),experimental_sem(1:maxlag), 'ok', 'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',2,'LineWidth',0.1);
    h.CapSize = 0.1;
    if  strcmp(geneFile , 'H2B_withTags.txt') == 1
        lineProps.col= {[1 0.6 0]};
        lineProps.width = 1;
        A = mseb(lags, simulation_autocorrelation,simulation_sem_autocorrelation,lineProps,0);
    end
    if  strcmp(geneFile , 'Bactin_withTags.txt') == 1
        lineProps.col={[0 .6 1]};
        lineProps.width = 1;
        A = mseb(lags, simulation_autocorrelation,simulation_sem_autocorrelation,lineProps,0);
    end
    if  strcmp(geneFile , 'KDM5B_withTags.txt') == 1
        lineProps.col = {[0.4 .0 1]};
        lineProps.width = 1;
        A = mseb(lags, simulation_autocorrelation,simulation_sem_autocorrelation,lineProps,0);
    end
    lgd=legend([h, A.mainLine],'Data','Model');
    set(lgd,'FontSize',5);
    lgd.Location='northeast';
    grid off;
    box on
    set(gca,'linewidth',0.5)
    xlabel('\tau (sec)','FontSize',12);
    ylabel('G(\tau)/G(0)','FontSize',12);
    xlim([0 250])
    ylim([-0.5 1.2])
    xticks([0 50 100 150 200 250])
    xticklabels({'0','50','100','150','200','250'})
    set (gca ,'FontSize',8, 'FontName', 'Arial');
    nameplot = horzcat('AutoCorr_',geneFile(1:end-13));
    print('-dpng','-r300',nameplot)
    movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
end

%% Evaluating Objective function
dataPoints = experimental_lags;
noDataPoints = length(dataPoints);
simulation_autocorrelation = simulation_autocorrelation(1:noDataPoints);
exp_var = diag(covariance_Autocorrelation)';
exp_var = exp_var(1:noDataPoints);
exp_var (exp_var==0)=1e-3;
fitAC =  (nTraces/2) * sum( ((experimental_autocorrelation(2:end) - simulation_autocorrelation(2:end)).^2)  ./ exp_var(2:end));

%  fitAC = fitAC/noDataPoints
%% New implementation
% load('kk.mat')
% N = 150
% matrix_Autocorrelation = matrix_Autocorrelation(:,2:N);
% SIGMA = cov(matrix_Autocorrelation);
% [U,S] = svd(SIGMA);
% Sred = S(1:Ntr-1,1:Ntr-1);
% Ured = U(:,1:Ntr-1);
% V=U';
% Vred = V(1:Ntr-1,:);
% SIGMAredinv = Ured*Sred^-1*Vred;
% dx = experimental_autocorrelation(2:N)- simulation_autocorrelation(2:N);
% true_AC = (Ntr/2)* (dx*SIGMAredinv*dx')
% fitAC =  (nTraces/2) * sum( ((experimental_autocorrelation(2:N) - simulation_autocorrelation(2:N)).^2)  ./ exp_var(2:N))

% save kk.mat
%  pause

end
