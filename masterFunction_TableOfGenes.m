clear all; close all;clc;
cd Source_Code

%% Code to creates an excel file with the genes used in this study.
%% This database contains information for 2647 genes from the PantherDB.

%% saving all genes
load ('Lib_Genes_nonRepeated.mat');
numElementsLibrary = length (Lib_unq);
for i =1:numElementsLibrary
    DB_Genes(i).Accession_gene = Lib_unq(i).Accession_gene;
    DB_Genes(i).geneLength = Lib_unq(i).geneLength;
    DB_Genes(i).geneLengthWithTag = Lib_unq(i).geneLength + 318;    
    DB_Genes(i).Sequence = Lib_unq(i).Sequence;
end
writetable(struct2table(DB_Genes), 'DB_Genes.xlsx')
cd ..