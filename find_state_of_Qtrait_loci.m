function [state_of_Qtrait_loci]=find_state_of_Qtrait_loci(loci) 

state_of_Qtrait_loci=cell(2,2);

%%find deleterious loci IN THE K QUANTITATIVE LOCI
idx_del = find(loci==1);
state_of_Qtrait_loci{1,2}= idx_del; % indices of deleterious loci, from 1 to K
state_of_Qtrait_loci{1,1}= size(idx_del,2); % number of deleterious loci

%%find benign loci
idx_ben = find(loci==0);
state_of_Qtrait_loci{2,2}= idx_ben;
state_of_Qtrait_loci{2,1}= size(idx_ben,2);








