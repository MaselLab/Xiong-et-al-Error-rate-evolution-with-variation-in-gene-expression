function [fixation_flux_of_drho,Max_Pfix]=P_drho(Delta,MutaRate_of_rho,POP_size,rho,s,Sigma_f,N_of_Qtrait_loci,trait_distance,state_of_Qtrait_loci,Current_exp_of_ldel,Tot_exp_frq,Range_drho,GridNum,Pben,Pdel)
global Grid_of_drho Alpha Beta PDF_of_drho 

Grid_of_drho_rescaled = 10.^Grid_of_drho * rho;
idx_of_possible_mutation=find(Grid_of_drho_rescaled<=1); %% make sure rho mutates to be greater than 1 is not possible
N_of_possible_mutation=idx_of_possible_mutation(end); 
Grid_of_drho_rescaled = Grid_of_drho_rescaled-rho;
Interval = 2*Range_drho/GridNum;

% this is the ratios of Wmisfold*Wproof_reading after and before mutations
s_of_drho_piece1=max(0,1-s/Tot_exp_frq*(Current_exp_of_ldel.*rho.*10.^Grid_of_drho(1:N_of_possible_mutation)+Pdel/(Pben+Pdel)*...
                (Tot_exp_frq-Current_exp_of_ldel)*rho^2*10.^(2*Grid_of_drho(1:N_of_possible_mutation))))/max(0,1-s/Tot_exp_frq*...
                (rho*Current_exp_of_ldel + Pdel/(Pben+Pdel)*(Tot_exp_frq-Current_exp_of_ldel)*rho^2)).*...
                (1-Delta.*log10(rho))./(1-Delta.*(log10(rho)+Grid_of_drho(1:N_of_possible_mutation)));    

if state_of_Qtrait_loci{2,1}==0  % if all quantative loci are deleterious, then changes in trait is determined by alpha and drho
    difference_by_drho = -sum(Alpha,1);
else % otherwise, beta needs to be considered
    operator_lben = zeros(N_of_Qtrait_loci,1);
    operator_lben(state_of_Qtrait_loci{2,2})=1;
    operator_lben = diag(operator_lben);
    difference_by_drho = sum(-Alpha+operator_lben*(Alpha+Beta),1);
end

s_of_drho = s_of_drho_piece1.*exp(-sum((Grid_of_drho_rescaled(1:N_of_possible_mutation)*difference_by_drho).^2 +...
            2.*Grid_of_drho_rescaled(1:N_of_possible_mutation)*(difference_by_drho.*trait_distance),2)/(2*Sigma_f^2))-1;

Pfix_of_drho = expm1(-s_of_drho)./expm1(-POP_size*s_of_drho); 

Pfix_of_drho(s_of_drho==0) = 1/POP_size;

fixation_flux_of_drho = MutaRate_of_rho*Interval*sum(Pfix_of_drho.*PDF_of_drho(1:N_of_possible_mutation));

Max_Pfix = max(Pfix_of_drho);