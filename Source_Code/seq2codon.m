function codons = seq2codon(seq)
% function converts a character array to cell array with each position
% containing codons
L = length(seq)/3; % codon length
if ceil(L) ~= floor(L)
    error('Invalid sequence length.')
end
codons = cell(1,L);
k = 1;
for i = 1:L
    codons{i} = seq(k:k+2);
    k = k+3;
end
end