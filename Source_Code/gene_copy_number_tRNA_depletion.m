% Gene Copy Number another version.
% 
GCN = struct(...
'ATA',  7.5, 'ATT', 16.0, 'ATC', 20.8,... % I Isoleucine
'CTA',  7.2, 'TTA',  7.7, 'TTG', 12.9,... % L Leucine
'CTT', 13.2, 'CTC', 19.6, 'CTG', 39.6,... % L Leucine
'GTA',  7.1, 'GTT', 11.0, 'GTC', 14.5,... % V Valine 
'GTG', 28.1,...                           % V Valine
'TTT', 17.6, 'TTC', 20.3,...              % F Phenylalanine
'ATG', 22.0,...                           % M START CODON (Methionine)
'TGT', 10.6, 'TGC', 12.6,...              % C Cysteine
'GCG', 7.4, 'GCA', 15.8, 'GCT', 18.4,...  % A Alanine
'GCC', 27.7,...                           % A Alanine
'GGT', 10.8, 'GGA', 16.5, 'GGG', 16.5,... % G Glycine
'GGC', 22.2,...                           % G Glycine
'CCG',  6.9, 'CCA', 16.9, 'CCT', 17.5,... % P Proline
'CCC', 19.8,...                           % P Proline
'ACG', 6.1, 'ACT', 13.1, 'ACA', 15.1,...  % T Threonine
'ACC', 18.9,...                           % T Threonine
'TCG',  4.4, 'AGT', 12.1, 'TCA', 12.2,... % S Serine 
'TCT', 15.2, 'TCC', 17.7, 'AGC', 19.5,... % S Serine
'TAT', 12.2, 'TAC', 15.3,...              % Y Tyrosine
'TGG', 13.2,...                           % W Tryptophan
'CAA', 12.3, 'CAG', 34.2,...              % Q Glutamine
'AAT', 17.0, 'AAC', 19.1,...              % N Asparganine
'CAT', 10.9, 'CAC', 15.1,...              % H Histidine
'GAA', 29.0, 'GAG', 39.6,...              % E Glutamic acid
'GAT', 21.8, 'GAC', 25.1,...              % D Aspartic acid
'AAA', 24.4, 'AAG', 31.9,...              % K Lysine
'CGT',  4.5, 'CGA',  6.2, 'CGC', 10.4,... % R Arginine
'CGG', 11.4, 'AGG', 12.0, 'AGA', 12.2,... % R Arginine
'TAG',  0.8, 'TAA',  1.0, 'TGA',  1.6);   % STOP codons 

GCN2 = struct(...
'ATA',  7.5, 'ATT', 16.0, 'ATC', 20.8,... % I Isoleucine
'CTA',  7.2, 'TTA',  7.7, 'TTG', 12.9,... % L Leucine
'CTT', 13.2, 'CTC', 19.6, 'CTG', 39.6,... % L Leucine
'GTA',  7.1, 'GTT', 11.0, 'GTC', 14.5,... % V Valine 
'GTG', 28.1,...                           % V Valine
'TTT', 17.6, 'TTC', 20.3,...              % F Phenylalanine
'ATG', 22.0,...                           % M START CODON (Methionine)
'TGT', 10.6, 'TGC', 12.6,...              % C Cysteine
'GCG', 7.4, 'GCA', 15.8, 'GCT', 18.4,...  % A Alanine
'GCC', 27.7,...                           % A Alanine
'GGT', 10.8, 'GGA', 16.5, 'GGG', 16.5,... % G Glycine
'GGC', 22.2,...                           % G Glycine
'CCG',  6.9, 'CCA', 16.9, 'CCT', 17.5,... % P Proline
'CCC', 19.8,...                           % P Proline
'ACG', 6.1, 'ACT', 13.1, 'ACA', 15.1,...  % T Threonine
'ACC', 18.9,...                           % T Threonine
'TCG',  4.4, 'AGT', 12.1, 'TCA', 12.2,... % S Serine 
'TCT', 15.2, 'TCC', 17.7, 'AGC', 19.5,... % S Serine
'TAT', 12.2, 'TAC', 15.3,...              % Y Tyrosine
'TGG', 13.2,...                           % W Tryptophan
'CAA', 12.3, 'CAG', 34.2,...              % Q Glutamine
'AAT', 17.0, 'AAC', 19.1,...              % N Asparganine
'CAT', 10.9, 'CAC', 15.1,...              % H Histidine
'GAA', 29.0, 'GAG', 39.6,...              % E Glutamic acid
'GAT', 21.8, 'GAC', 25.1,...              % D Aspartic acid
'AAA', 24.4, 'AAG', 31.9,...              % K Lysine
'CGT',  4.5, 'CGA',  6.2, 'CGC', 10.4,... % R Arginine
'CGG', 11.4, 'AGG', 12.0, 'AGA', 12.2); % R Arginine



minmaxAA = struct(...
'I', [7.5 20.8], 'L', [7.2 39.6], 'V', [7.1 28.1], 'F', [17.6 20.3],...
'M', [22 22], 'C', [10.6 12.6], 'A', [7.4 27.7], 'G', [10.8 22.2],...
'P', [6.9 19.8], 'T', [6.1 18.9], 'S', [4.4 19.5], 'Y', [12.2 15.3],...
'W', [13.2 13.2], 'Q', [12.3 34.2], 'N', [17 19.1], 'H', [10.9 15.1],...
'E', [29 39.6], 'D', [21.8 25.1], 'K', [24.4 31.9], 'R', [4.5 12.2]);

fields = fieldnames(GCN);
% si = 1; fi = 1;
for i = 1:numel(fields)
    if GCN.(fields{i}) < 10
        sGCN.(fields{i}) = GCN.(fields{i});
    elseif GCN.(fields{i}) > 30
        fGCN.(fields{i}) = GCN.(fields{i});
    end
end

codonAA = struct(...
'ATA', 'I', 'ATT', 'I', 'ATC', 'I',...  % I Isoleucine
'CTA', 'L', 'TTA', 'L', 'TTG', 'L',...  % L Leucine
'CTT', 'L', 'CTC', 'L', 'CTG', 'L',...  % L Leucine
'GTA', 'V', 'GTT', 'V', 'GTC', 'V',...  % V Valine 
'GTG', 'V',...                          % V Valine
'TTT', 'F', 'TTC', 'F',...              % F Phenylalanine
'ATG', 'M',...                          % M START CODON (Methionine)
'TGT', 'C', 'TGC', 'C',...              % C Cysteine
'GCG', 'A', 'GCA', 'A', 'GCT', 'A',...  % A Alanine
'GCC', 'A',...                          % A Alanine
'GGT', 'G', 'GGA', 'G', 'GGG', 'G',...  % G Glycine
'GGC', 'G',...                          % G Glycine
'CCG', 'P', 'CCA', 'P', 'CCT', 'P',...  % P Proline
'CCC', 'P',...                          % P Proline
'ACG', 'T', 'ACT', 'T', 'ACA', 'T',...  % T Threonine
'ACC', 'T',...                          % T Threonine
'TCG', 'S', 'AGT', 'S', 'TCA', 'S',...  % S Serine 
'TCT', 'S', 'TCC', 'S', 'AGC', 'S',...  % S Serine
'TAT', 'Y', 'TAC', 'Y',...              % Y Tyrosine
'TGG', 'W',...                          % W Tryptophan
'CAA', 'Q', 'CAG', 'Q',...              % Q Glutamine
'AAT', 'N', 'AAC', 'N',...              % N Asparganine
'CAT', 'H', 'CAC', 'H',...              % H Histidine
'GAA', 'E', 'GAG', 'E',...              % E Glutamic acid
'GAT', 'D', 'GAC', 'D',...              % D Aspartic acid
'AAA', 'K', 'AAG', 'K',...              % K Lysine
'CGT', 'R', 'CGA', 'R', 'CGC', 'R',...  % R Arginine
'CGG', 'R', 'AGG', 'R', 'AGA', 'R',...  % R Arginine
'TAG', 'STOP', 'TAA', 'STOP', 'TGA', 'STOP');   % STOP codons 



AA.I = {'ATA', 'ATT', 'ATC'};
AA.L = {'CTA','TTA','TTG','CTT','CTC','CTG'};
AA.V = {'GTA','GTT','GTC','GTG'};
AA.F = {'TTT','TTC'};
AA.M = {'ATG'};
AA.C = {'TGT','TGC'};
AA.A = {'GCG','GCA','GCT','GCC'};
AA.G = {'GGT','GGA','GGG','GGC'};
AA.P = {'CCG','CCA','CCT','CCC'};
AA.T = {'ACG','ACT','ACA','ACC'};
AA.S = {'TCG','AGT','TCA','TCT','TCC','AGC'};
AA.Y = {'TAT','TAC'};
AA.W = {'TGG'};
AA.Q = {'CAA','CAG'};
AA.N = {'AAT','AAC'};
AA.H = {'CAT','CAC'};
AA.E = {'GAA','GAG'};
AA.D = {'GAT','GAC'};
AA.K = {'AAA','AAG'};
AA.R = {'CGT','CGA','CGC','CGG','AGG','AGA'};

% GCN = struct(...
% 'UUU', 17.6,  'UCU', 15.2,  'UAU', 12.2,  'UGU', 10.6,...
% 'UUC', 20.3,  'UCC', 17.7,  'UAC', 15.3,  'UGC', 12.6,...
% 'UUA',  7.7,  'UCA', 12.2,  'UAA',  1.0,  'UGA',  1.6,...
% 'UUG', 12.9,  'UCG',  4.4,  'UAG',  0.8,  'UGG', 13.2,...
% 'CUU', 13.2,  'CCU', 17.5,  'CAU', 10.9,  'CGU',  4.5,...
% 'CUC', 19.6,  'CCC', 19.8,  'CAC', 15.1,  'CGC', 10.4,...
% 'CUA',  7.2,  'CCA', 16.9,  'CAA', 12.3,  'CGA',  6.2,...
% 'CUG', 39.6,  'CCG',  6.9,  'CAG', 34.2,  'CGG', 11.4,...
% 'AUU', 16.0,  'ACU', 13.1,  'AAU', 17.0,  'AGU', 12.1,...
% 'AUC', 20.8,  'ACC', 18.9,  'AAC', 19.1,  'AGC', 19.5,...
% 'AUA',  7.5,  'ACA', 15.1,  'AAA', 24.4,  'AGA', 12.2,...
% 'AUG', 22.0,  'ACG',  6.1,  'AAG', 31.9,  'AGG', 12.0,...
% 'GUU', 11.0,  'GCU', 18.4,  'GAU', 21.8,  'GGU', 10.8,...
% 'GUC', 14.5,  'GCC', 27.7,  'GAC', 25.1,  'GGC', 22.2,...
% 'GUA',  7.1,  'GCA', 15.8,  'GAA', 29.0,  'GGA', 16.5,...
% 'GUG', 28.1,  'GCG',  7.4,  'GAG', 39.6,  'GGG', 16.5);

