function plotTimeCourses (time, intensityVector,nExpRepetitions,geneFile,folderName,samplingRateInSeconds)


for i =1: size (intensityVector,1)
 intensityVector(i,:) = intensityVector(i,:) ./ max(intensityVector(i,:));
%intensityVector(i,:) = (intensityVector(i,:) - mean(intensityVector(i,:)))' / std(intensityVector(i,:)) ;

end

time = downsample(time,samplingRateInSeconds);
intensityVector = downsample (intensityVector', samplingRateInSeconds)';

selectedSpot=6;

figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
% fig1.PaperPosition = [0, 0, 2.3, 1.1/2];
fig1.PaperPosition = [0, 0, 2.3, 0.6];

for i =1:nExpRepetitions
    plot(time,intensityVector(i,:),'-','Color',[0.7,0.7,0.7],'LineWidth',0.2);
    hold on
end
if  strcmp(geneFile , 'H2B_withTags.txt') == 1
    plot(time,intensityVector(selectedSpot,:),'Color',[1 0.6 0],'LineWidth',1);
 %   ylim([0 1])
end
if  strcmp(geneFile , 'Bactin_withTags.txt') == 1
    plot(time,intensityVector(selectedSpot,:),'Color',[0 .6 1],'LineWidth',1);
 %   ylim([0 1])
% xlim([0 800])
end
if  strcmp(geneFile , 'KDM5B_withTags.txt') == 1
    plot(time,intensityVector(selectedSpot,:),'Color',[0.4 .0 1],'LineWidth',1);
 %   ylim([0 1])
%
end

 yticks([0 0.5 1])
yticklabels({'0','0.5','1'})

 xlim([0 250])
 ylim([0 1])
 
grid off;
box on
set(gca,'linewidth',1)
xlabel('time (sec)','FontSize',10);
ylabel('I(t) (au)','FontSize',10);
set (gca ,'FontSize',6, 'FontName', 'Arial');
nameplot = horzcat('TC_',geneFile(1:end-13));
print('-dpng','-r300',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');


close all;
