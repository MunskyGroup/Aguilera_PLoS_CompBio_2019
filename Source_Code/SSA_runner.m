function [RibosomePositions,L] = SSA_runner(sequenceFile,t,nR,el_rate,in_rate,opt,speed)
%function [k_e,Lags,mean_AC,err_AC,t_out,cell_I,tRNA_cn,RibosomePositions,L] = SSA_runner(codons,t,nR,el_rate,in_rate,opt,speed)
%% 1.- Obtaining gene sequence from GeneBank
% [~, ~,probePos,tagcdn] = sequenceAnalyzer(seq);
% tag_info;
% codons = [tagcdn, codons];
% L = size(codons,2);
%% 2.- Defining Parameters for simulation
k_e_Mean = el_rate;
k_i = in_rate;

%% 3.- Running Stochastic Model
if (~exist('opt','var')) && (~exist('speed','var'))
    [L,RibosomePositions] = control_Fn(sequenceFile, k_e_Mean, k_i, nR, t);
elseif (exist('opt','var')) && (~exist('speed','var'))
    [L, RibosomePositions] = control_Fn(sequenceFile, k_e_Mean, k_i, nR, t, opt);
elseif (exist('opt','var')) && (exist('speed','var'))
    [L,RibosomePositions] = control_Fn(sequenceFile, k_e_Mean, k_i, nR, t, opt, speed);
else
    error('Improper input.')
end
end