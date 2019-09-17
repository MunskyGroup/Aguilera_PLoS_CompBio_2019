function plot_Collisions (probability_Collision,p_space1,p_space2,nameplot)

%% Plotting
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 3.4, 2.5];
surf(p_space1,p_space2,probability_Collision);
view(2)
colormap('jet');
hcb=colorbar;
lim = caxis;
caxis([0 max(max(probability_Collision))]);
xlim([min(p_space1),max(p_space1)]);
ylim([min(p_space2),max(p_space2)]);
% ylabel(hcb, '# Collisions / Ribosome','FontSize',12);
ylabel(hcb, 'No Collisions','FontSize',12);

ylabel('k_{e} (aa/sec)','FontSize',10);
xlabel('k_{i} (1/sec)','FontSize',10);
% nameplot = horzcat('parameterScan_Collisions_median');
set (gca ,'FontSize',8, 'FontName', 'Arial');
print('-dpng','-r600',nameplot)

end
