function [L, RibosomePositions] = control_Fn(geneFile,k_eMean,k_i,nR,t_sim,opt,speed)
% codons = cell array of all codons
% k_eMean = assumed average elongation rate
% k_i = assumed initiation rate
% probePos = vector containing codon positions of the probes
% nR = number of repetitions for simulations
% t_sim = time to which simulation will run
% opt = codon of interest (COI). Example: 'CAG'. optional
% speed = allows for optimization or deoptimization or tRNA knock down/in
% 'slow' = changes to slowest synonymous codon for codon of interest
% 'fast' = changes to fastest synonymous codon for codon of interest
% 'decrease' = 0.05x the gene copy number for the codon of interest
% 'increase' = 20x the gene copy number for the codon of interest
% # = multiply tRNA copy number by number input for 'speed'

% Outputs:
% L = gene length in codons
% time = output times from SSA
% cell_I = cell intensity vector from SSA
% ribDens = ribosome density
% tRNA_cn = tRNA 


[generated_geneSequence,codons] = sequenceAnalyzer_tRNA_Depletion(geneFile);

gene_copy_number_tRNA_depletion;
% M = 1 selects for the slowest codon for an AA.
% M = 2 selects for the fastest codon for an AA.
if (~exist('speed','var'))
    % if a speed is not provided, set the default to slow, which is M = 1.
    M = 1;
elseif strcmp(speed,'slow') == 1
    M = 1;
elseif strcmp(speed,'fast') == 1
    M = 2;
end

if (~exist('opt','var'))
    % if no options exist, default to running the natural sequence
    ver = 0;
elseif exist('opt','var')
    % if opt exists, a codon of interest has been specified. 
    ver = 1;
    cdn = opt;
    AAi = codonAA.(cdn); % determines which aa specified codon codes for
    minmax = minmaxAA.(AAi); % the min/max gene copy numbers for COI.
    % the first value of minmax is the min, second is the max
    if exist('M','var') && minmax(M) == GCN.(cdn)
        warning('Copy number of selected codon is already at a minimum or maximum.')
    end
end

try
if strcmp(speed,'decrease') == 1
    % set multiplying factor f to 0.05
    ver = 2;
    f = 0.05; 
elseif strcmp(speed,'increase') == 1
    % set multiplying factor f to 20
    ver = 2;
    f = 20;
elseif isnumeric(speed)
    ver = 2;
    f = speed;
elseif exist('speed','var')
    error('Invalid input for speed option.')
end
catch
    
end


%% Calculating the size of the nt and aa sequences.
L = length(codons);

%% Elongation constant.
k_elongation=  zeros (1,L );
tRNA_cn = zeros (1,L );

% Ver 0 supplies natural unmodified gene copy numbers for [tRNA]
if ver == 0
 %    disp('Ver = 0')
    for i = 1 : L
        field = codons{i};
        tRNA_cn(i) = GCN.(field) ;
    end   
end

% Ver 1 changes the 
if ver == 1
  %   disp('Ver = 1')
    for i = 1 : L
        field = codons{i};
        if strcmp(field,cdn) == 1
            tRNA_cn(i) = minmax(M);
        else
            tRNA_cn(i) = GCN.(field);
        end   
    end
end

% Ver 2 either increases or decrease [tRNA] for the specified codon
if ver == 2
     %disp('Ver = 2')
    for i = 1 : L
        field = codons{i};
        if strcmp(field,cdn) == 1
            tRNA_cn(i) = f*GCN.(field);
        else
            tRNA_cn(i) = GCN.(field);
        end
    end   
end


% calculating the mean values in the structure
%mean_tRNACopyNumber = mean (mean(reshape(struct2array(GCN),numel(fieldnames(GCN)),[]),2));
mean_tRNACopyNumber = mean (mean(reshape(struct2array(GCN2),numel(fieldnames(GCN2)),[]),2));


for i = 1 : L
    k_elongation(i)= (tRNA_cn (i)/ mean_tRNACopyNumber)* k_eMean;
end

all_Constants = [k_i,k_elongation,10];

nonConsiderTime = 0;

%% creating the model
% nPoints= t_sim;
% probePositionVector = zeros (1, L);
% for i =1: length(probePos)
%     probePositionVector(1, probePos(i):end) = 1*i;
% end
% time = linspace(0,t_sim,t_sim);


TimeVectorFixedSize = linspace (0,nonConsiderTime+t_sim,nonConsiderTime+t_sim);

parfor k = 1 : nR
    % Running the methods
%    [RibosomePositions{1,k}] = solve_direct_optimized(nonConsiderTime,nonConsiderTime+t_sim, all_Constants);

%   RibosomePositions{1,k} = SSA_longNames_mex(all_Constants,TimeVectorFixedSize,0,0);
  RibosomePositions{1,k}  = SSA_longNames_mex(all_Constants,TimeVectorFixedSize,0,0,0);

  
  

%     s_out_Short {1,k}= s_out(nonConsiderTime:end,:);
%     No_ribosomes(k,:) = sum ( s_out(nonConsiderTime:end,:));
end
% No_ribosomes = mean (No_ribosomes);
% ribDens = No_ribosomes/nPoints;
% rib_per_RNA = mean (No_ribosomes);
% 
% %% Saving the intensity vectors in a cell array.
% intensityVector = zeros (nR,t_sim );
% I  = zeros (1,t_sim );
% for k = 1 : nR
%     for tp =1: t_sim
%         temp_simulationOutPut = s_out_Short{1,k}(tp,1:L);
%         I(tp)= sum (probePositionVector .* temp_simulationOutPut );
%     end
%     intensityVector (k,:)= I;
% end
% cell_I {1,1} = intensityVector;

end
