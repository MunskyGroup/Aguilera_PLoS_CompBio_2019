function plottingIntensityDistribution(hisData,geneFile, folderName)
hisData = hisData./10;
%% Plotting Distributions
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 2.5, 1.6];

hold on
if  strcmp(geneFile , 'H2B_withTags.txt') == 1
    colorPlot =[1 0.6 0];
end

if strcmp(geneFile ,'Bactin_withTags.txt') == 1
    colorPlot =[0 .6 1];
end

if strcmp(geneFile ,'KDM5B_withTags.txt') == 1
    colorPlot =[0.4 .0 1];
end

colorPlot_OP= [1 0 0];
colorPlot_DeOP= [0 0 0];
maXInt =15;

for i =1:3
    edg= 0:1:maXInt;
    hisData_sim_norm =histcounts ( hisData(i,:),edg,'Normalization','probability');
    switch i
        case 1
            h(i) = stairs (hisData_sim_norm,'-', 'Color',colorPlot, 'LineWidth',2);
        case 2
            h(i) = stairs (hisData_sim_norm,':*', 'Color',colorPlot_OP, 'LineWidth',2);
        case 3
            h(i) = stairs (hisData_sim_norm,'-.s', 'Color',colorPlot_DeOP, 'LineWidth',1);
    end
    hold on
end
lgd=legend([h(1),h(2),h(3)],{'Natural','Common', 'Rare'});
set(lgd,'FontSize',8);
xlim([1 , 15])
box on
set(gca,'linewidth',1)
grid off
set (gca ,'FontSize',8);
set(gca, 'FontName', 'Arial');
xlabel('Intensity (ump)','FontName', 'Arial','FontSize',8);
ylabel('Probability','FontName', 'Arial','FontSize',8);
nameplot = horzcat('Dis_Int_',geneFile(1:end-13));
print('-dpng','-r300',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');

end