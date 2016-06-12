function [fixation_flux_of_dalpha,fixation_flux_of_dbeta, max_Pfix_of_dalpha, max_Pfix_of_dbeta]=P_alpha_beta_biased(MutRate_of_alpha,MutRate_of_beta,POP_size,rho,Sigma_f,N_of_Qtrait_loci,trait_distance,state_of_Qtrait_loci,N_of_Qtrait,Range_dalpha,Range_dbeta,GridNum,a)
global PDF_of_dalpha PDF_of_dbeta Grid_of_dalpha Grid_of_dbeta Alpha Beta

        %NOTE this code is not suitable for multiple Qtrait. Because one
        %quantitative locus can contribute to multiple traits, a mutation
        %to it should affect multiple trait-specific alpha and beta. Here, 
        %only one alpha or one beta is mutated at a time, which can only
        %happen if there is only one quantitative trait. For this reason,
        %N_of_Qtrait is essentially fixed at one.
        
I1 = ones(1,GridNum);
I2 = ones(N_of_Qtrait*N_of_Qtrait_loci,1);

%% dalpha
I4=ones(N_of_Qtrait_loci,N_of_Qtrait);
I4(state_of_Qtrait_loci{1,2},:)=0;
I5=reshape(I4,N_of_Qtrait_loci*N_of_Qtrait,1);

trait_dist_reshaped_1 = reshape(ones(N_of_Qtrait_loci,1)*trait_distance,N_of_Qtrait*N_of_Qtrait_loci,1);

Alpha_reshaped = reshape(Alpha,N_of_Qtrait*N_of_Qtrait_loci,1);

diff_by_dalpha = (1-rho)*(Alpha_reshaped*I1./(-a)+I2*Grid_of_dalpha)+rho*((I5.*Alpha_reshaped)*I1/(-a)+I5*Grid_of_dalpha); % alpha goes row-wise, dalpha goes column-wise.

s_of_dalpha = exp(-(diff_by_dalpha.^2+2*diff_by_dalpha.*(trait_dist_reshaped_1*I1))./(2*Sigma_f^2))-1; % this is the selection coefficient

Pfix_of_dalpha = expm1(-s_of_dalpha)./expm1(-s_of_dalpha.*POP_size);

Pfix_of_dalpha(s_of_dalpha==0) = 1/POP_size;

Pfix_of_dalpha = sum(Pfix_of_dalpha,1); % this sums pfix over loci and traits, yielding the overall pfix of each dalpha

max_Pfix_of_dalpha = max(Pfix_of_dalpha); 

Interval_Alpha = 2*Range_dalpha/GridNum;

fixation_flux_of_dalpha = MutRate_of_alpha*Interval_Alpha*sum(PDF_of_dalpha.*Pfix_of_dalpha); 


%% dbeta
Interval_Beta = 2*Range_dbeta/GridNum;

if state_of_Qtrait_loci{2,1}~= 0 %% calc the contribution from lben only when there are lben in the K set
    
    I3 = ones(state_of_Qtrait_loci{2,1}*N_of_Qtrait,1);
        
    trait_dist_reshaped_2 = reshape(ones(state_of_Qtrait_loci{2,1},1)*trait_distance,N_of_Qtrait*state_of_Qtrait_loci{2,1},1);

    Beta_reshaped = reshape(Beta(state_of_Qtrait_loci{2,2},:),state_of_Qtrait_loci{2,1}*N_of_Qtrait,1);
    
    s_of_dbeta_lben = exp(-((rho.*(Beta_reshaped*I1./(-a) + I3*Grid_of_dbeta)).^2 + 2.*rho.*(Beta_reshaped*I1./(-a) + I3*Grid_of_dbeta).*(trait_dist_reshaped_2*I1))./(2*Sigma_f^2))-1;   

    Pfix_of_dbeta_lben = expm1(-s_of_dbeta_lben)./expm1(-s_of_dbeta_lben.*POP_size);

    Pfix_of_dbeta_lben(s_of_dbeta_lben==0) = 1/POP_size;  
        
    Pfix_of_dbeta_lben = sum(Pfix_of_dbeta_lben,1);
    
    if state_of_Qtrait_loci{1,1}~=0 %% calc the contribution from ldel only when there are ldel in the K set
        Pfix_of_dbeta_ldel = state_of_Qtrait_loci{1,1}.*N_of_Qtrait./POP_size.*ones(1,GridNum); % for a deleterious locus, mutations to beta have no effect on fitness, so pfix=1/POP_size

        Pfix_of_dbeta = Pfix_of_dbeta_lben+Pfix_of_dbeta_ldel;
        
        max_Pfix_of_dbeta = max(Pfix_of_dbeta);

        fixation_flux_of_dbeta = MutRate_of_beta*Interval_Beta*sum(PDF_of_dbeta.*Pfix_of_dbeta); 
    else
        Pfix_of_dbeta = Pfix_of_dbeta_lben;

        max_Pfix_of_dbeta = max(Pfix_of_dbeta_lben);
        
        fixation_flux_of_dbeta = MutRate_of_beta*Interval_Beta*sum(Pfix_of_dbeta.*PDF_of_dbeta);

    end

else
    Pfix_of_dbeta = N_of_Qtrait_loci*N_of_Qtrait./POP_size;
    
    max_Pfix_of_dbeta = Pfix_of_dbeta;
    
    fixation_flux_of_dbeta = MutRate_of_beta*Interval_Beta*sum(PDF_of_dbeta.*Pfix_of_dbeta);
end

