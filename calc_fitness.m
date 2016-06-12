function [fitness,current_Qtrait,Wmisfold,trait_distance]=calc_fitness(rho,s,Sigma_f,Delta,state_of_Qtrait_loci,N_of_Qtrait_loci,N_of_Qtrait,Tot_exp_frq,Opt_trait,Exp_frq, Pben, Pdel)
global Loci Alpha Beta

Current_exp_of_ldel = Loci*Exp_frq;

Wmisfold = max(0,1-s/Tot_exp_frq*(Current_exp_of_ldel*rho+(Tot_exp_frq-Current_exp_of_ldel)*rho^2*Pdel/(Pben+Pdel))); 

operator = zeros(N_of_Qtrait_loci,1);
    
if state_of_Qtrait_loci{2,1}~= 0
    operator(state_of_Qtrait_loci{2,2})=1;
end

current_x_matrix = (1-rho)*Alpha + rho*diag(operator)*(Alpha+Beta);
current_Qtrait= sum(current_x_matrix,1);
trait_distance = current_Qtrait-Opt_trait;
Wtrait = (1/(2*pi)^0.5/Sigma_f)^N_of_Qtrait*exp(-sum(trait_distance.^2)/(2*Sigma_f^2)); % the coefficient (1/(2*pi)^0.5/Sigma_f) is left out in the article, but it will not affect our result.
fitness = Wmisfold*Wtrait/(1-Delta*log10(rho));
