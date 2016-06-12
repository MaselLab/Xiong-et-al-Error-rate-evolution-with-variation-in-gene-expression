function [id_in_set_K,N_of_ldel, state_of_Qtrait_loci]=muta_coopt(state_of_Qtrait_loci, SD_of_NewBeta,rho,trait_distance,...
                                                                  N_of_ldel, fixation_flux_ben, fixation_flux_ben_per_locus, max_Pfix_of_NewBeta_ben, ...
                                                                  fixation_flux_del, fixation_flux_del_per_locus,Sigma_f,POP_size,Ltot,N_of_Qtrait_loci)

global Alpha Beta Loci idx_of_lben idx_of_ldel

% step 1: is the new cryptic sequence benign or deleterious?
P_ben=fixation_flux_ben/(fixation_flux_ben+fixation_flux_del);
ref=rand(1);
if ref<=P_ben
   tag=1;
else
   tag=0;
end

switch tag
    case 1 % benign
        weight=fixation_flux_ben_per_locus./sum(fixation_flux_ben_per_locus);
        id_in_set_K = datasample(state_of_Qtrait_loci{2,2},1,'Weight',weight);  % step 2: which locus   
        ceiling=1.1*max_Pfix_of_NewBeta_ben(state_of_Qtrait_loci{2,2}==id_in_set_K);
        % step 3: what value for the new beta
        while 1
            ref=rand(1)*ceiling;
            NewBeta=normrnd(0,SD_of_NewBeta);
            difference_by_beta = (1-rho)*Beta(id_in_set_K)+rho*NewBeta; 
            s=exp(-(difference_by_beta^2+2*trait_distance*difference_by_beta)/(2*Sigma_f^2))-1;
            if s~=0
                Pfix=expm1(-s)./expm1(-s*POP_size);
            else
                Pfix=1/POP_size;
            end
            if ref<=Pfix
                break;
            end           
        end        
    case 0 % deleterious
        weight=fixation_flux_del_per_locus./sum(fixation_flux_del_per_locus);
        id_in_set_K = datasample(state_of_Qtrait_loci{2,2},1,'Weight',weight);
        NewBeta=normrnd(0,SD_of_NewBeta);      
end

%% now, apply this cooption
Alpha(id_in_set_K)=Alpha(id_in_set_K)+Beta(id_in_set_K); 
Beta(id_in_set_K)=NewBeta; 

if tag==0 % if the cryptic sequence is deleterious, also update loci and state_of_Qtrait_loci
    id_in_loci = id_in_set_K+Ltot-N_of_Qtrait_loci;
    Loci(id_in_loci)=1;
    idx_of_lben(idx_of_lben==id_in_loci)=[];
    idx_of_ldel = [idx_of_ldel,id_in_loci];
    N_of_ldel = N_of_ldel+1;
    state_of_Qtrait_loci{2,1} = state_of_Qtrait_loci{2,1}-1;
    state_of_Qtrait_loci{2,2}(state_of_Qtrait_loci{2,2}==id_in_set_K) = [];
    state_of_Qtrait_loci{1,1} = state_of_Qtrait_loci{1,1}+1;
    state_of_Qtrait_loci{1,2} = [state_of_Qtrait_loci{1,2},id_in_set_K];
end



