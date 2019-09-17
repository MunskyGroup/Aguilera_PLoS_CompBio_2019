function plotGenomeDistributions(ke_hat,geneLength, namePlot,minXval,maxXval,folderName)

ke_hat(isnan(ke_hat)==1)=0;

numElementsLibrary = length(ke_hat);

for i =1:numElementsLibrary
    if geneLength(1,i)>500 && geneLength(1,i)<= 1000
        category(i) = 2; % medium
    elseif geneLength(1,i) >1000
        category(i) = 3; % long
    elseif geneLength(1,i)<=500
        category(i) = 1; % short
    end
end
edg=linspace(5,20,30);
[N_s,edges_s] = histcounts (ke_hat(category==1),edg,'Normalization','probability');
[N_m,edges_m] = histcounts (ke_hat(category==2),edg,'Normalization','probability');
[N_l,edges_l] = histcounts (ke_hat(category==3),edg,'Normalization','probability');
edges_s = edges_s(2:end) - (edges_s(2)-edges_s(1))/2;
edges_m = edges_m(2:end) - (edges_m(2)-edges_m(1))/2;
edges_l = edges_l(2:end) - (edges_l(2)-edges_l(1))/2;

figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 2.3*1.1, 1.2*1.1];
box on
hold on
l1 = plot (edges_s,N_s,'-m ','LineWidth',2,'Color',[1 0.6 0]);
p1 = area(edges_s,N_s,'FaceColor',[1 0.6 0]);
p1.FaceAlpha = 0.2;
% childp1=get(p1,'Children');
% set(childp1,'FaceAlpha',0.01)
%set(childp1,'FaceColor',[1 0.6 0]) % 0 .6 1 % 0.4 .0 1

l2 = plot (edges_m,N_m,'-b','LineWidth',1,'Color',[0 .6 1]);
p2 = area(edges_m,N_m,'FaceColor',[0 .6 1]);
p2.FaceAlpha = 0.2;
%childp2=get(p2,'Children');
%set(childp2,'FaceAlpha',0.01)

l3 = plot (edges_l,N_l,'-r','LineWidth',1,'Color',[0.4 .0 1]);
p3 = area(edges_l,N_l,'FaceColor',[0.4 .0 1]);
p3.FaceAlpha = 0.2;
% childp3=get(p3,'Children');
% set(childp3,'FaceAlpha',0.01)

p4 =plot ([10,10],[0, 1.2],'--k','LineWidth',1);
%lgd =legend([l1,l2,l3,p4],'Short','Medium','Long','ke =10');
%set(lgd,'FontSize',5);
% set(lgd,'Location','northwest');
%set(lgd,'Location','northeast');
%title('Association Time')

xlabel('ke (aa/sec)','FontSize',12);
ylabel('Probability','FontSize',12);

xlim([minXval,maxXval])
ylim([0,0.55])

% xticks([0.2 0.4])
% xticklabels({'0.2','0.4'})
% yticks([5 10 15])
% yticklabels({'5','10','15'})

print('-dpng','-r300',namePlot)
movefile(horzcat(namePlot, '.png'),horzcat(folderName),'f');


end