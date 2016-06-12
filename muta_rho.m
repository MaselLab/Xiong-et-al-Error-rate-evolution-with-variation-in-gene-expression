function [rho]=muta_rho(Delta,POP_size,rho,s,Sigma_f,Current_exp_of_ldel,Tot_exp_frq,state_of_Qtrait_loci,trait_distance,N_of_Qtrait_loci,Max_Pfix,Pben,Pdel,bias_mut_rho)
global Alpha Beta

if state_of_Qtrait_loci{2,1}==0
    difference_by_drho = -sum(Alpha,1);
else
    operator_lben = zeros(N_of_Qtrait_loci,1);
    operator_lben(state_of_Qtrait_loci{2,2})=1;
    operator_lben = diag(operator_lben);
    difference_by_drho = sum(-Alpha+operator_lben*(Alpha+Beta),1);
end

while 1        
    temp=10; % no magic about this value. see the next couple lines.    
    while temp>1 % make sure rho <=1 after mutation
        drho = normrnd(0, 0.2)+bias_mut_rho;
        temp = 10^drho*rho;
    end    
    
    drho_rescaled = (10^drho-1)*rho;   

    s_of_drho_piece1 = max(0,1-s/Tot_exp_frq*(rho*10^drho*Current_exp_of_ldel + Pdel/(Pdel+Pben)*(Tot_exp_frq-Current_exp_of_ldel)*(rho*10^drho)^2))...
                       /max(0,1-s/Tot_exp_frq*(rho*Current_exp_of_ldel + Pdel/(Pdel+Pben)*(Tot_exp_frq-Current_exp_of_ldel)*rho^2)).*...
                       (1-Delta.*log10(rho))./(1-Delta.*(log10(rho)+drho));

    s_of_drho = s_of_drho_piece1.*exp(-sum((drho_rescaled*difference_by_drho).^2 + 2.*drho_rescaled*(difference_by_drho.*trait_distance),2)/(2*Sigma_f^2))-1;

    if s_of_drho == 0
        Pfix_of_drho = 1/POP_size;
    else
        Pfix_of_drho = expm1(-s_of_drho)./expm1(-POP_size*s_of_drho);      
    end
    
    % below is the rejection method
    reference = rand(1)*Max_Pfix*1.1;    
    if Pfix_of_drho >= reference;
        break;
    end
end
    
rho = rho*10^drho;