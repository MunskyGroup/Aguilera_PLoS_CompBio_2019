function gene_sequence_visualizer2(gene_seq, opt,design,folderName)
close all
gene_copy_number_tRNA_depletion;
tag_info; % generates tagcdn


sCDN = fieldnames(sGCN); % slow codon names
fCDN = fieldnames(fGCN); % fast codon names
if (~exist('opt','var'))
    opt = sCDN;
    col = [1 0 0]; % color for slow codons
elseif strcmp(opt,'fast') == 1
    opt = fCDN;
    col = [0 0 1]; % color for fast codons
elseif strcmp(opt,'slow') == 1
    opt = sCDN;
    col = [1 0 0]; 
elseif exist('opt','var') == 1
    opt = opt;
    col = [0 0 0]; % color for selected codon
end
grey = [0.9, 0.9,0.9];

% seq should be a character string
seq = upper(gene_seq);
if iscell(gene_seq) == 0
    geneLength = length(seq)/3;
    
    if geneLength ~= floor(geneLength) || geneLength ~= ceil(geneLength)
        error('Invalid gene length.')
    end
    codons = zeros(0,geneLength);
    p = 1;
    for i = 1:geneLength
        codons{i} = seq(p:p+2);
        p = p+3;
    end
end

if iscell(gene_seq) == 1
    codons = gene_seq;
    geneLength = length(codons);
end


%% Plotting
figure('visible', 'on');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, 3, 0.25];
hold on

% h1=rectangle('Position', [length(tagcdn) 0 (geneLength-length(tagcdn))  1],'Facecolor','w','EdgeColor','k');
% h1.LineWidth = 0.1;

%k=0;

h2=rectangle('Position', [0 0 length(tagcdn) 1], 'Facecolor',grey,'Edgecolor','k');
h2.LineWidth = 0.1;

plot ([1,1],[0,1],'k-','LineWidth',0.25)    

for i = 1:geneLength
    if any(strcmp(opt,codons{i})) == 1
     %   k = k+1;
     % plot ([length(tagcdn)+i-1,length(tagcdn)+i-1],[0,1],'k-','LineWidth',0.25)    
      plot ([i-1,i-1],[0,1],'k-','LineWidth',0.25)  
    else
%         fc = [1 ,1 ,1];
    end
end
xlim([1,geneLength])
box on
set(gca,'linewidth',0.25)
set(gcf, 'Position', [50, 300, 1200, 100])
 set(gca, 'YTick',[], 'XTick',[])
%xlim([0 geneLength+length(tagcdn)])
%disp(k)
nameplot = horzcat('Design_',design);
print('-dpng','-r600',nameplot)
movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');

%k

end