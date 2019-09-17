function [codonUsage, codonSensitivity, CAI_codons, CAI] = sequenceCodonAnalysis(geneSequence)
codonUsage = zeros(1,20);
%% Calculating the size of the nt and aa sequences.
geneSequence = upper(geneSequence);
codons = geneSequence;
geneLength = length(codons)/3;
proteinFrame = nt2aa(geneSequence);

%% Codon frequency.
codonUsage(1,1) = length(strfind(proteinFrame,'A'));
codonUsage(1,2) = length(strfind(proteinFrame,'R'));
codonUsage(1,3) = length(strfind(proteinFrame,'N'));
codonUsage(1,4) = length(strfind(proteinFrame,'D'));
codonUsage(1,5) = length(strfind(proteinFrame,'C'));

codonUsage(1,6) = length(strfind(proteinFrame,'Q'));
codonUsage(1,7) = length(strfind(proteinFrame,'E'));
codonUsage(1,8) = length(strfind(proteinFrame,'G'));
codonUsage(1,9) = length(strfind(proteinFrame,'H'));
codonUsage(1,10) = length(strfind(proteinFrame,'I'));

codonUsage(1,11) = length(strfind(proteinFrame,'L'));
codonUsage(1,12) = length(strfind(proteinFrame,'K'));
codonUsage(1,13) = length(strfind(proteinFrame,'M'));
codonUsage(1,14) = length(strfind(proteinFrame,'F'));
codonUsage(1,15) = length(strfind(proteinFrame,'P'));

codonUsage(1,16) = length(strfind(proteinFrame,'S'));
codonUsage(1,17) = length(strfind(proteinFrame,'T'));
codonUsage(1,18) = length(strfind(proteinFrame,'W'));
codonUsage(1,19) = length(strfind(proteinFrame,'Y'));
codonUsage(1,20) = length(strfind(proteinFrame,'V'));

%% Normalizing Codon usage with respect to the geneLength
codonUsageNormalized = codonUsage./ geneLength;

%% Calling the section of the code that reads the codon usage and ratioFast_Slow
genomicCopyNumber;

%% Deffining the sensitivity function
codonSensitivity = codonUsageNormalized .*  ratio_Sentitivity_Fast_Slow;
codonSensitivity = round (codonSensitivity,2);

%% Calculating the CAI
counter = 1;
for i =1 : geneLength-1
    separated_codons(i,:) = codons(counter : counter + 2);
    counter = counter + 3;
end



for i = 1 : geneLength-1
    field = separated_codons(i,:);
    CAI_codons (i) = strGenCopy.(field) / strGenCopy_Fast.(field) ;
end

CAI = geomean (CAI_codons);
end
