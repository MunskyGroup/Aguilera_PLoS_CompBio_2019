clear all; close all; clc;
cd Source_Code

%% This code is inteded to simulate single molecule stochastic dynamics.
%% The code reproduces figure 3 center panel.
%% All results are stored in Results_Methods.

folderName = horzcat('Results_Methods'); if exist (folderName, 'dir') ~= 7; mkdir(folderName);end
%% Compiling a single library.
load ('Lib_short.mat');
numElementsLibrary = length (Lib_unq_short);
ki=0.03;
asuumed_ke =10;
meanHT_Delay =10;
%% Prealocating memory
vector_NoSpots = [10, 30, 60, 100];
N_Rep_Spots = length (vector_NoSpots);

%% Running the SSA for the real Ke.
ke_SR1 = zeros (numElementsLibrary,N_Rep_Spots);ke_SR3=ke_SR1;ke_SR10=ke_SR1;ke_SR20=ke_SR1;
geneLength = zeros (1,numElementsLibrary);
for i =1: numElementsLibrary
    %       HT = i
    [generated_geneSequence, typeOfTag,tagPositions,parametersModel ,geneLength(i)] = sequenceAnalyzer(Lib_unq_short(i).Sequence,ki,asuumed_ke);
    [ke_SR1(i,:),ke_SR3(i,:),ke_SR10(i,:),ke_SR20(i,:)] = function_calculate_HT_ke_REP_2 (parametersModel,tagPositions,vector_NoSpots,meanHT_Delay);
end

%% calculating codon usage for each gene.
codonUsage = zeros (1,numElementsLibrary);
for  i =1: numElementsLibrary
    tagSeq='ATGGACTACAAGGACGACGACGACAAAGGTGACTACAAAGATGATGACGATAAAGGCGACTATAAGGACGATGACGACAAGGGCGGAAACTCACTGATCAAGGAAAACATGCGGATGAAGGTGGTGATGGAGGGCTCCGTGAATGGTCACCAGTTCAAGTGCACCGGAGAGGGAGAGGGAAACCCGTACATGGGAACTCAGACCATGCGCATTAAGGTCATCGAAGGAGGTCCGCTGCCGTTCGCTTTCGATATCCTGGCCACTTCGTTCGGAGGAGGGTCGCGCACGTTCATCAAGTACCCGAAGGGAATCCCGGACTTCTTTAAGCAGTCATTCCCGGAAGGATTCACTTGGGAACGGGTGACCCGGTATGAAGATGGAGGTGTGGTGACTGTCATGCAAGATACTTCGCTGGAGGATGGGTGCCTCGTGTACCACGTCCAAGTCCGCGGAGTGAATTTCCCGTCCAACGGACCAGTGATGCAGAAAAAGACGAAGGGTTGGGAACCTAATACTGAAATGATGTACCCCGCAGACGGAGGGCTGAGGGGCTACACCCACATGGCGCTGAAGGTCGACGGAGGAGATTACAAGGATGACGACGATAAGCAACAAGATTACAAAGACGATGATGACAAGGGCCAGCAGGGCGACTACAAGGACGACGACGACAAGCAGCAGGACTACAAAGATGACGATGATAAAGGAGGAGGACATCTGTCCTGTTCGTTCGTGACCACCTACAGATCAAAGAAAACCGTGGGAAACATCAAGATGCCGGGCATTCATGCCGTCGACCACCGCCTGGAGCGGCTCGAAGAATCAGACAATGAGATGTTCGTCGTGCAAAGAGAACATGCCGTGGCCAAGTTCGCGGGACTGGGAGGCGGTGGAGGCGATTACAAAGACGATGATGACAAGGGTGACTATAAAGACGACGATGACAAAGGGGATTACAAGGATGATGATGATAAGGGAGGCGGTGGATCAGGTGGAGGAGGTTCACTGCAG';
    geneSequence = strcat(tagSeq,Lib_unq_short(i).Sequence);
    codonUsage (i) = sequenceUsage(geneSequence);
end

%% calculate AC
ke_SR1c = zeros (numElementsLibrary,N_Rep_Spots); ke_SR3c=ke_SR1c;ke_SR10c=ke_SR1c;ke_SR20c=ke_SR1c;
for  i =1: numElementsLibrary
    ke_SR1c(i,:) =  ke_SR1(i,:)  .*(1/ geneLength(i)) .*codonUsage(i);
    ke_SR3c(i,:) =  ke_SR3(i,:)  .*(1/ geneLength(i)) .*codonUsage(i);
    ke_SR10c(i,:) = ke_SR10(i,:) .*(1/ geneLength(i)) .*codonUsage(i);
    ke_SR20c(i,:) = ke_SR20(i,:) .*(1/ geneLength(i)) .*codonUsage(i);
end

%% Calculating RMSE.
[RMSE_long(1,:),RMSE_medium(1,:),RMSE_short(1,:)] = function_RMSE (ke_SR1c, N_Rep_Spots,numElementsLibrary,geneLength,asuumed_ke); % correlations sampling rate =1
[RMSE_long(2,:),RMSE_medium(2,:),RMSE_short(2,:)] = function_RMSE (ke_SR3c, N_Rep_Spots,numElementsLibrary,geneLength,asuumed_ke);  % correlations sampling rate =3
[RMSE_long(3,:),RMSE_medium(3,:),RMSE_short(3,:)] = function_RMSE (ke_SR10c, N_Rep_Spots,numElementsLibrary,geneLength,asuumed_ke);  % correlations sampling rate =5
[RMSE_long(4,:),RMSE_medium(4,:),RMSE_short(4,:)] = function_RMSE (ke_SR20c, N_Rep_Spots,numElementsLibrary,geneLength,asuumed_ke);  % correlations sampling rate =10

%% Saving data
save HT_ke_REP1.mat ke_SR1 ke_SR3 ke_SR10 ke_SR20 ke_SR1c ke_SR3c ke_SR10c ke_SR20c RMSE_long RMSE_medium RMSE_short numElementsLibrary geneLength asuumed_ke

%% Plotting results
% load HT_ke_REP
namePlot ='HT_r';
plot_ElongationMethods (folderName, namePlot, RMSE_long, RMSE_medium, RMSE_short)

%% Plotting Distributions
minXval =5; maxXval =20;
nameplot ='dist_HT';
plotGenomeDistributions(ke_SR3c(:,1),geneLength, nameplot,minXval,maxXval,folderName)

cd ..