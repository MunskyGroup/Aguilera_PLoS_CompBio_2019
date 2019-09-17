close all; clc; clear all;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics with tRNA depletion stress.
%% The code reproduces figure S7.
%% All results are stored in Results_tRNA_Depletion.

folderName = horzcat('Results_tRNA_Depletion'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end
geneFile1 = 'H2B_withTags.txt'; inRate(1)=0.066;
geneFile2 = 'Bactin_withTags.txt'; inRate(2)=0.05;
geneFile3 = 'KDM5B_withTags.txt';  inRate(3)=0.022;

%% Running the SSA
N = 20; % number of tRNA concentrations tested
trna_f=exp(linspace(log(0.001),log(10),N)); % logspace between tRNA concentrations
tSim = 100000;
nR = 8;
elRate = 10.6;
%inRate = 0.033;
ke1 = zeros(1,N);
ke2 = ke1; ke3 = ke1;
% associaiton time = AT
AT1 = zeros(1,N);
AT2 = AT1; AT3 = AT1;
CDN_depleted ='CTC';

nRepetitions =3;
%% plotting sequence visualization
[~,codons] = sequenceAnalyzer_tRNA_Depletion(geneFile1);
gene_sequence_visualizer2(codons,CDN_depleted,'H2B',folderName);
[~,codons] = sequenceAnalyzer_tRNA_Depletion(geneFile2);
gene_sequence_visualizer2(codons,CDN_depleted,'Bact',folderName);
[~,codons] = sequenceAnalyzer_tRNA_Depletion(geneFile3);
gene_sequence_visualizer2(codons,CDN_depleted,'KDM5B',folderName);


%% Simulating the model
for k =1:nRepetitions
    for i = 1:N
        [RibosomePositions,L] = SSA_runner('H2B_withTags.txt',tSim,nR,elRate,inRate(1),CDN_depleted,trna_f(i));
        [ke1(k,i),AT1(k,i)] = realElongationRates(RibosomePositions{1,1},L);

        [RibosomePositions,L] = SSA_runner('Bactin_withTags.txt',tSim,nR,elRate,inRate(2),CDN_depleted,trna_f(i));
        [ke2(k,i),AT2(k,i)] = realElongationRates(RibosomePositions{1,1},L);
        
        [RibosomePositions,L] = SSA_runner('KDM5B_withTags.txt',tSim,nR,elRate,inRate(3),CDN_depleted,trna_f(i));
        [ke3(k,i),AT3(k,i)] = realElongationRates(RibosomePositions{1,1},L);
    end
end

%save simData.mat ke1 ke2 ke3 trna_f AT1 AT2 AT3

colors = [241, 90, 90;
    240, 196, 25;
    78, 186, 111;
    45, 149, 191;
    149, 91, 165]/255;

%% Plotting
close all
%figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
% fig1.PaperPosition = [0, 0, 2.2, 2.5];
fig1.PaperPosition = [0, 0, 2.8, 2];
hold on
lineProps.col= {[1 0.6 0]};
lineProps.width = 1.5;
d1 = mseb(trna_f,mean(ke1),std(ke1),lineProps,1);
lineProps.col={[0 .6 1]};
lineProps.width = 1.5;
d2 = mseb(trna_f,mean(ke2),std(ke2),lineProps,1);
lineProps.col = {[0.4 .0 1]};
lineProps.width = 1.5;
d3 = mseb(trna_f,mean(ke3),std(ke3),lineProps,1);
p1=plot(ones(1,10),linspace(0,12,10),'-r','LineWidth',1);
p2=plot(ones(1,10)*1e-2,linspace(0,12,10),'--r','LineWidth',1);
xlabel('[tRNA_{CTC}]','FontSize',12)
ylabel('k_{e} (aa/sec)','FontSize',12)
lgd=legend([d1.mainLine,d2.mainLine,d3.mainLine,p1,p2],'H2B','b-act', 'KDM5B','Naural [tRNA_{CTC}]','99% depletion');
set(lgd,'FontSize',5);
lgd.Location='southeast';
box on
set(gca,'XScale','log')
set(gca,'linewidth',1)
xlim([min(trna_f) max(trna_f)]);
ylim([0 12]);
set (gca ,'FontSize',8, 'FontName', 'Arial');
nameplot = 'tRNA_concentration';
print('-dpng','-r600',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');

%% Plotting Association time
close all
%figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 2.8, 2];
% fig1.PaperPosition = [0, 0, 2.2, 2.5];
hold on
lineProps.col= {[1 0.6 0]};
lineProps.width = 2;
d1 = mseb(trna_f,mean(AT1),std(AT1),lineProps,1);
lineProps.col={[0 .6 1]};
lineProps.width = 2;
d2 = mseb(trna_f,mean(AT2),std(AT2),lineProps,1);
lineProps.col = {[0.4 .0 1]};
lineProps.width = 2;
d3 = mseb(trna_f,mean(AT3),std(AT3),lineProps,1);
p1=plot(ones(1,10),linspace(10,1e5,10),'-r','LineWidth',1);
p2=plot(ones(1,10)*1e-2,linspace(10,1e5,10),'--r','LineWidth',1);
xlabel('[tRNA_{CTC}]','FontSize',12)
ylabel('Association Time (sec)','FontSize',12)
lgd=legend([d1.mainLine,d2.mainLine,d3.mainLine,p1,p2],'H2B','b-act', 'KDM5B','Naural [tRNA_{CTC}]','99% depletion');
set(lgd,'FontSize',5);
lgd.Location='northeast';
box on
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'linewidth',1)
xlim([min(trna_f) max(trna_f)]);
%ylim([0 11]);
set (gca ,'FontSize',8, 'FontName', 'Arial');
nameplot = 'tRNA_concentration_AT';
print('-dpng','-r600',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');

cd ..
