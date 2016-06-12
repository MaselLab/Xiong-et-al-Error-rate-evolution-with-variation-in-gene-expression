function [fixation_flux_of_coopt, fixation_flux_ben, fixation_flux_ben_per_locus, max_Pfix_of_NewBeta_ben, fixation_flux_del, fixation_flux_del_per_locus]=...
    P_coopt(rho,Sigma_f,state_of_Qtrait_loci,trait_distance,POP_size,Pdel,Pben,Current_exp_of_ldel,Current_Wm, Range_of_NewBeta, Ltot, N_of_Qtrait_loci,GridNum,N_of_Qtrait,s,Tot_exp_frq)
global Alpha Beta Grid_of_NewBeta PDF_of_NewBeta Loci Trans_exp_frq

interval=2*Range_of_NewBeta/GridNum;

% the idea is to draw a coopt event with three steps in muta_coopt. 
% step 1: is the new cryptic seq ben or del
% step 2: which locus 
% step 3: what value of the new beta 
% Therefore, fixation fluxes of the above types must be calculated here

if state_of_Qtrait_loci{2,1}==0;    % if all quantitative loci are deleterious, then no cooption
    fixation_flux_of_coopt=0;   
    fixation_flux_del=0;
    fixation_flux_ben=0;
    fixation_flux_ben_per_locus=0;
    fixation_flux_del_per_locus=0;
    max_Pfix_of_NewBeta_ben=0;   
else
    I1 = ones(1,GridNum);    
    I3 = ones(state_of_Qtrait_loci{2,1}*N_of_Qtrait,1);       
    trait_dist_reshaped_2 = reshape(ones(state_of_Qtrait_loci{2,1},1)*trait_distance,N_of_Qtrait*state_of_Qtrait_loci{2,1},1); 

    Beta_reshaped = reshape(Beta(state_of_Qtrait_loci{2,2},:),state_of_Qtrait_loci{2,1}*N_of_Qtrait,1); 
    
    %% if co-option creats a benign locus
    delta_trait_ben = (1-rho).*Beta_reshaped*I1+rho.*I3*Grid_of_NewBeta;  % if the new locus is benign, the change in trait is (1-rho)*beta+rho*new_beta
    
    ratio_of_Wtrait_ben = exp(-(delta_trait_ben.^2+2.*delta_trait_ben.*(trait_dist_reshaped_2*I1))./(2*Sigma_f^2)); % new beta are placed column-wise, beta are placed row-wise
    
    s_of_coopt_piece1 = ratio_of_Wtrait_ben-1;  
    
    Pfix_of_coopt_piece1 = expm1(-s_of_coopt_piece1)./expm1(-s_of_coopt_piece1.*POP_size);
    
    Pfix_of_coopt_piece1(s_of_coopt_piece1==0) = 1/POP_size;    
    
    fixation_fluxes_of_coopt_piece1 = interval.*Pfix_of_coopt_piece1.*(I3*PDF_of_NewBeta);
       
    fixation_flux_ben = sum(sum(fixation_fluxes_of_coopt_piece1)).*Pben./(Pben+Pdel);   % fixation flux of all cooptions that produce benign cryptic new seq. To be used in step 1
    
    fixation_flux_ben_per_locus = sum(fixation_fluxes_of_coopt_piece1,2);  % locus-wise fixation flux of all cooptions that produce benign cryptic new seq. To be used in step 2. Note that the calculation at this step is correct only when N_of_Qtrait=1
    
    max_Pfix_of_NewBeta_ben = max(Pfix_of_coopt_piece1,[],2); % these are the locus-wise max pfix of new beta 
    
    %% if the new cryptic sequence is deleterious, then things are complicate  
    Alpha_reshaped = reshape(Alpha(state_of_Qtrait_loci{2,2},:),state_of_Qtrait_loci{2,1}*N_of_Qtrait,1);
    
    delta_trait_del = (1-2*rho).*Beta_reshaped-rho.*Alpha_reshaped; % the change in trait is (1-2*rho)*beta-rho*alpha
    
    ratio_of_Wtrait_del = exp(-(delta_trait_del.^2+2.*delta_trait_del.*trait_dist_reshaped_2)./(2*Sigma_f^2)); % this is just a vector
    
    ratio_of_Wm = (max(0,1-s/Tot_exp_frq*(Current_exp_of_ldel*rho+(Tot_exp_frq-Current_exp_of_ldel)*rho^2*Pdel/(Pdel+Pben))-s/Tot_exp_frq*...
                  (rho-rho^2*Pdel/(Pdel+Pben))*((Loci(state_of_Qtrait_loci{2,2}+Ltot-N_of_Qtrait_loci)+1)...
                  .*Trans_exp_frq(state_of_Qtrait_loci{2,2}+Ltot-N_of_Qtrait_loci))))/Current_Wm;
                 
    s_of_coopt_piece2 = ratio_of_Wtrait_del.*reshape(ratio_of_Wm*ones(1,N_of_Qtrait),state_of_Qtrait_loci{2,1}*N_of_Qtrait,1)-1;   %% dim=(state_of_Qtrait_loci{2,1}*N_of_Qtrait,1), compatiable with N_of_Qtrait>1
    
    Pfix_of_coopt_piece2 = expm1(-s_of_coopt_piece2)./expm1(-s_of_coopt_piece2.*POP_size);
    
    Pfix_of_coopt_piece2(s_of_coopt_piece2==0) = 1/POP_size;
    
    fixation_fluxes_of_coopt_piece2 = interval.*Pfix_of_coopt_piece2*PDF_of_NewBeta;
    
    fixation_flux_del = sum(sum(fixation_fluxes_of_coopt_piece2)).*Pdel./(Pben+Pdel);   % to be used in step 1
    
    fixation_flux_del_per_locus = sum(fixation_fluxes_of_coopt_piece2,2);  % to be used in step 2     
    
    %% combine the cases of new cryptic sequence being benign and being deleterious
    fixation_flux_of_coopt=fixation_flux_ben+fixation_flux_del;    
end

