function plotCollisionDistributions_ki (collisions,pre_namePlot,folderName)
close all

%% plotting
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 2.5, 2];
fig_pos = fig1.PaperPosition;
fig1.PaperSize = [fig_pos(3) fig_pos(4)];
max_x_val = max(collisions);
N = histcounts(collisions, 'BinLimits', [0 max_x_val], 'BinMethod', 'integers', 'Normalization', 'probability');
xvalues = linspace (0,max_x_val,max_x_val+1)-0.5;
% h= stairs(xvalues,N);
h= histogram (collisions,40);
% h(1).Marker = 'o';
% h(1).MarkerSize = 4;
% h(1).MarkerFaceColor = 'm';
h(1).FaceColor = [0.4,0.4,0.4];
h(1).LineWidth = 2;
h(1).Normalization = 'probability';
%xlim([-0.5,max_x_val]);
xlim([-0.5, 200]);
%ylim([0,1.1]);
box on
set(gca,'linewidth',1)
set (gca ,'FontSize',10);
set(gca, 'FontName', 'Arial');
xlabel('No Collisions','FontSize',10,'FontName', 'Arial');
ylabel('Probability','FontName', 'Arial','FontSize',10);
nameplot = horzcat('Prob_',pre_namePlot);
print('-dpng','-r300',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
end