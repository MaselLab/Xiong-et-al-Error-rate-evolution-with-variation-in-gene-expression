function muta_alpha_beta_biased(rho,POP_size,Sigma_f,Sigma_m,trait_distance,N_of_Qtrait_loci,N_of_Qtrait,Max_Pfix,a,Ltot,type)

global Alpha Beta Loci

I1 = ones(N_of_Qtrait_loci,N_of_Qtrait);
trait_dist_matrix = ones(N_of_Qtrait_loci,1)*trait_distance;
I2 = (1-Loci(Ltot-N_of_Qtrait_loci+1:Ltot)')*ones(1,N_of_Qtrait);

switch type
    case 1 %% mut alpha
        %the strategy here is to first pick a mutation that happens equally 
        %for each alpha, i.e equally possible for each loci and each trait,
        % and check the overall probability of fixation of this mutation.
        %If it fixes, draw a locus based on the locus specific probability 
        %of fixation. 
        
        %The computation of the ratio of Wtrait after and before mutations is made unnecessarily complicate (line26-line39).        
        %Note the exponent in ratio is[(1-rho)*(dalpha-alpha/a)+I2*rho*(dalpha-alpha/a)]^2+2*[(1-rho)*(dalpha-alpha/a)+I2*rho*(dalpha-alpha/a)]*trait_distance (ignoring the -1/2*sigma_f^2)
        %Expanding and rearranging the terms, we have
        %[(1-rho)*(-alpha/a)+I2*rho*(-alpha/a)]^2+2*[(1-rho)*(-alpha/a)+I2*rho*(-alpha/a)]*trait_distance+[(1-rho)*dalpha+I2*rho*dalpha]^2+2*[(1-rho)*dalpha+I2*rho*dalpha]*[(1-rho)*(-alpha/a)+I2*rho*(-alpha/a)+trait_distance]
        %Calling [(1-rho)*(-alpha/a)+I2*rho*(-alpha/a)] Alpha_revised and (1-rho)+I2*rho I3, we have
        %Alpha_revised^2+2*Alpha_revised*trait_distance+(I3.*dalpha)^2+2*I3.*dalpha.*(Alpha_revised+trait_distance)
        %The first two terms is processed first and named ratio_of_Wtrait_piece1.
        %Similarly for mutations to beta (line61-71).        
        
       Alpha_revised = (1-rho)./(-a).*Alpha + rho/(-a)*I2.*Alpha;
       
       ratio_of_Wtrait_piece1 = exp(-(Alpha_revised.^2 + 2*Alpha_revised.*trait_dist_matrix)./(2*Sigma_f^2));
       
       Alpha_plus_trait_dist = Alpha_revised + trait_dist_matrix;
       
       I3=(1-rho)*I1+rho*I2;
              
       while 1
            reference = rand(1)*Max_Pfix*1.1;

            dalpha = normrnd(0, Sigma_m/N_of_Qtrait_loci);

            s_of_dalpha = ratio_of_Wtrait_piece1.*exp(-((I3.*dalpha).^2+2.*I3.*dalpha.*Alpha_plus_trait_dist)./(2*Sigma_f^2))-1;           

            Pfix_of_dalpha_matrix = expm1(-s_of_dalpha)./expm1(-s_of_dalpha.*POP_size);

            Pfix_of_dalpha_matrix(s_of_dalpha==0) = 1/POP_size;
            
            Pfix_of_dalpha = sum(sum(Pfix_of_dalpha_matrix));
  
            if Pfix_of_dalpha >= reference;
                break;
            end
        end
        
        weight = Pfix_of_dalpha_matrix./Pfix_of_dalpha;

        weight = reshape(weight,1,N_of_Qtrait_loci*N_of_Qtrait);
        
        idx = datasample(1:N_of_Qtrait_loci*N_of_Qtrait,1,'Weight',weight);
                    
        Alpha(idx)=Alpha(idx).*(a-1)/a + dalpha;
           
    case 2
        Beta_revised = rho*Beta.*I2./(-a);
        Beta_plus_trait_dist = Beta_revised + trait_dist_matrix;
        ratio_of_Wtrait_piece1 = exp(-(Beta_revised.^2 + 2*Beta_revised.*trait_dist_matrix)./(2*Sigma_f^2));
        
        while 1
         
            reference = rand(1)*Max_Pfix*1.1;

            dbeta =normrnd(0, Sigma_m/N_of_Qtrait_loci); 
            
            s_of_dbeta = ratio_of_Wtrait_piece1 .* exp(-(rho^2*(dbeta.*I2).^2 + 2*rho*dbeta.*I2.*Beta_plus_trait_dist)./(2*Sigma_f^2)) -1;        

            Pfix_of_dbeta_matrix = expm1(-s_of_dbeta)./expm1(-s_of_dbeta.*POP_size);

            Pfix_of_dbeta_matrix(s_of_dbeta==0) = 1/POP_size;
            
            Pfix_of_dbeta = sum(sum(Pfix_of_dbeta_matrix));

            if Pfix_of_dbeta>= reference
                break;
            end
        end
        
        weight = Pfix_of_dbeta_matrix./Pfix_of_dbeta;

        weight = reshape(weight,1,N_of_Qtrait_loci*N_of_Qtrait);
        
        idx = datasample(1:N_of_Qtrait_loci*N_of_Qtrait,1,'Weight',weight);
                    
        Beta(idx)=Beta(idx).*(a-1)/a + dbeta;
       
end

