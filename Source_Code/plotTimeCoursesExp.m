function plotTimeCoursesExp (fileName,maxlag,nExpRepetitions,geneFile,folderName,samplingRateInSeconds)


[~,sheet_name]=xlsfinfo(fileName);
for i=1:numel(sheet_name)
        pre_intensityVector =[];
    %% Extracting data and removing zeros
    data{i} = xlsread(fileName,sheet_name{i});
    temp1 = data{i}(:,2);
    StrData(1,i).trajectory = temp1;
    IntensityNormalized = temp1' / max(temp1);

    % IntensityNormalized = temp1';
    % IntensityNormalized = IntensityNormalized-mean(IntensityNormalized);
    % IntensityNormalized = IntensityNormalized/std(IntensityNormalized);
    if length(IntensityNormalized)>= maxlag
        pre_intensityVector(1,:) = IntensityNormalized(1:maxlag);
    else
        pre_intensityVector(1,:) = IntensityNormalized;
        pre_intensityVector(1,maxlag) = 0;
    end
        intensityVector(i,:) = pre_intensityVector;
end

time = linspace(1,maxlag,maxlag)*samplingRateInSeconds;
time(1)=0;
selectedSpot=8;

%% Figures
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
% fig1.PaperPosition = [0, 0, 2.3, 1.1/2];
fig1.PaperPosition = [0, 0, 2.3, 0.6];

for i =1:nExpRepetitions
    plot(time,intensityVector(i,:),'-','Color',[0.7,0.7,0.7],'LineWidth',0.2);
    hold on
    plot(time,intensityVector(selectedSpot,:),'Color',[0 0 0],'LineWidth',0.5);
end

 xlim([0 250])
 ylim([0 1]) 
 yticks([0 0.5 1])
yticklabels({'0','0.5','1'})
 
grid off;
box on
set(gca,'linewidth',1)
xlabel('time (sec)','FontSize',10);
ylabel('I(t) (au)','FontSize',10);
set (gca ,'FontSize',6, 'FontName', 'Arial');
nameplot = horzcat('TC_exp_',geneFile(1:end-13));
print('-dpng','-r300',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');


close all;
