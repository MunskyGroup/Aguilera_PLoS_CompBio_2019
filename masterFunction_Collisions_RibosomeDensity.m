clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces figure 5 top panel in the publication.
%% All results are stored in Results_Collisions_Sweep
%% For the Paper, b-act gene is selected.

folderName = horzcat('Results_Collisions_Sweep'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end

%% Deffing parameter values
nRepetitions = 8;
nonConsideredInitialSimulationTime = 2000;
totalSimulationTime = 2000+ nonConsideredInitialSimulationTime;
Resolution =30;
timePerturbationApplication= 0;
evaluatingInhibitor = 0;
evaluatingFRAP =0;
ke_vec = linspace(2,13,Resolution);
ki_vec = linspace(0.005,0.08,Resolution);
nPoints = totalSimulationTime - nonConsideredInitialSimulationTime;

%% loading experimental data
geneFile = 'Bactin_withTags.txt'; % 'H2B_withTags.txt'; 'KDM5B_withTags.txt';

%% Initiation Values
numberOfParameters =Resolution;
maximumValue1 = 0.08;
minimumValue1 =0.005;
p_space_ki = linspace( minimumValue1,maximumValue1 , Resolution );

%% Elongation Values
maximumValue2 = 2;
minimumValue2 =13;
p_space_ke = linspace( minimumValue2,maximumValue2 , Resolution );


%% Generating parameters for multiple ki and ke
parametersModel{numberOfParameters,numberOfParameters} =[];
for i=1: numberOfParameters
    for j=1: numberOfParameters
        [~, ~,~,parametersModel{i,j}] = sequenceAnalyzer(geneFile, p_space_ki (j), p_space_ke(i));
    end
end
geneLength = length(parametersModel{1,1})-2;

%% Calculating ribosome density
ribosomeDensity = zeros (numberOfParameters,numberOfParameters);
s_out{1,nRepetitions} = [];

for i=1: numberOfParameters
    for j=1: numberOfParameters
        %% Creating Model from Gene sequence
        parfor k = 1 : nRepetitions
            [~,s_out{k},~,~] = solve_direct_optimized(totalSimulationTime, parametersModel{i,j}, timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP);
            pre_ribosomeDensity(k,:) = sum (s_out{k}(nonConsideredInitialSimulationTime+1:totalSimulationTime,:),2)';
        end
        ribosomeDensity (i,j) = mean(mean(pre_ribosomeDensity));
    end
end

%% Calculating Number of bases between ribosomes
BasesPerRibosome = geneLength./ribosomeDensity;
[Xgrid1,Ygrid1] = meshgrid(p_space_ki,p_space_ke);

%% Saving files
fileToSave = horzcat('final_ribosomeDensity_',geneFile(1:end-13),'.mat');
save(fileToSave,'Xgrid1','Ygrid1','ribosomeDensity','geneFile','geneLength', 'BasesPerRibosome')

%% Plotting
figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 3.4, 2.5];
pcolor(Xgrid1,Ygrid1,BasesPerRibosome);
%contourf(Xgrid1,Ygrid1,BasesPerRibosome)
%shading interp
hcb=colorbar;
colormap('jet');
set(gca,'colorscale','log')
%colormap(flipud(colormap))
ylabel(hcb, 'Ribosome Distance','FontSize',12);
ylabel('k_{e} (aa/sec)','FontSize',10);
xlabel('k_{i} (1/sec)','FontSize',10);
nameplot = horzcat('CodonsPerRibosome_',geneFile(1:end-13));
set (gca ,'FontSize',8, 'FontName', 'Arial');
print('-dpng','-r600',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');

cd ..
