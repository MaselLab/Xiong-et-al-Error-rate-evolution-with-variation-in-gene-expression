function [Pfix_of_ldel,Pfix_of_lben]=P_ben_and_del(POP_size,rho,s,Sigma_f,N_of_Qtrait_loci,Ltot,trait_distance,state_of_Qtrait_loci,Current_Wmisfold,Current_exp_of_ldel,N_of_ldel,Tot_exp_frq,Pben,Pdel)
global Loci Trans_exp_frq idx_of_lben idx_of_ldel Beta Alpha
Pfix_of_ldel =0;
Pfix_of_lben =0;

%% ben to del
if N_of_ldel < Ltot  

    Wmisfold_of_new_ldel=max(0,1-s/Tot_exp_frq*(Current_exp_of_ldel*rho+(Tot_exp_frq-Current_exp_of_ldel)*rho^2*Pdel/(Pben+Pdel))-...
                   s/Tot_exp_frq*(rho-rho^2*Pdel/(Pben+Pdel))*((Loci(idx_of_lben)+1).*Trans_exp_frq(idx_of_lben)));

    if state_of_Qtrait_loci{2,1} ~= 0 %% if no lben in the K set, don't bother calculating the second piece      
        Ratio_of_Wtrait = zeros(1,Ltot); 
        
        Ratio_of_Wtrait(idx_of_lben) = 1;  
        
        difference_by_new_ldel = - rho.* diag(1-Loci(Ltot-N_of_Qtrait_loci+state_of_Qtrait_loci{2,2}))*(Alpha(state_of_Qtrait_loci{2,2},:)+Beta(state_of_Qtrait_loci{2,2},:));  
        
        Ratio_of_Wtrait(state_of_Qtrait_loci{2,2}+Ltot-N_of_Qtrait_loci)= exp(-sum(difference_by_new_ldel.^2 +...
                                                                       2*ones(state_of_Qtrait_loci{2,1},1)*trait_distance.*difference_by_new_ldel,2)/(2*Sigma_f^2));
        Ratio_of_Wtrait = (nonzeros(Ratio_of_Wtrait))';
        
        s_of_ldel = Wmisfold_of_new_ldel./Current_Wmisfold.*Ratio_of_Wtrait-1; %calc the overall selection coefficient of a new ldel
    else        
        s_of_ldel = Wmisfold_of_new_ldel./Current_Wmisfold -1;
    end     
  
    Pfix_of_ldel = expm1(-s_of_ldel)./expm1(-POP_size.*s_of_ldel);
    Pfix_of_ldel(s_of_ldel==0) = 1/POP_size;
end

%% del to ben
if N_of_ldel > 0
 
    Wmisfold_of_new_lben=max(0,1-s/Tot_exp_frq*(Current_exp_of_ldel*rho+(Tot_exp_frq-Current_exp_of_ldel)*rho^2*Pdel/(Pben+Pdel))...
                   +s/Tot_exp_frq*(rho-rho^2*Pdel/(Pben+Pdel))*(Loci(idx_of_ldel).*Trans_exp_frq(idx_of_ldel)));

    if state_of_Qtrait_loci{1,1} ~= 0 %% no ldel in K set, don't bother
        
        Ratio_of_Wtrait = zeros(Ltot,1);
        
        Ratio_of_Wtrait(idx_of_ldel) = 1;

        difference_by_new_lben = rho.* diag(Loci(Ltot-N_of_Qtrait_loci+state_of_Qtrait_loci{1,2}))*(Alpha(state_of_Qtrait_loci{1,2},:)+Beta(state_of_Qtrait_loci{1,2},:));

        Ratio_of_Wtrait(state_of_Qtrait_loci{1,2}+Ltot-N_of_Qtrait_loci)= exp(-sum(difference_by_new_lben.^2 + 2*ones(state_of_Qtrait_loci{1,1},1)*trait_distance.*difference_by_new_lben,2)/(2*Sigma_f^2));

        Ratio_of_Wtrait = (nonzeros(Ratio_of_Wtrait))';        
     
        s_of_lben = Wmisfold_of_new_lben/Current_Wmisfold.*Ratio_of_Wtrait-1;
    else
        
        s_of_lben = Wmisfold_of_new_lben./Current_Wmisfold -1;
        
    end

    Pfix_of_lben = expm1(-s_of_lben)./expm1(-POP_size.*s_of_lben);

    Pfix_of_lben(s_of_lben==0) = 1/POP_size;

end


