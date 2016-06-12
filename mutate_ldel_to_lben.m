function [N_of_ldel,state_of_Qtrait_loci,recalc_Wq] = mutate_ldel_to_lben(Ltot,N_of_Qtrait_loci, Pfix_of_lben, N_of_ldel, state_of_Qtrait_loci)
global Loci idx_of_ldel idx_of_lben

weight = Pfix_of_lben./sum(Pfix_of_lben);

idx = datasample(idx_of_ldel,1,'Weight',weight);

Loci(idx) = 0;

recalc_Wq = 1;

if idx < Ltot-N_of_Qtrait_loci+1   %% if the mutation is outside the K set, no need to recalculate Wtrait in the next round 
    
    recalc_Wq = 0;
    
end

idx_of_ldel(idx_of_ldel==idx)=[]; %% remove this index from the record of ldel

idx_of_lben = [idx_of_lben,idx]; %% add this index to the record of lben

N_of_ldel = N_of_ldel-1;

if idx > Ltot-N_of_Qtrait_loci

    state_of_Qtrait_loci{1,1}= state_of_Qtrait_loci{1,1}-1;

    state_of_Qtrait_loci{1,2}(state_of_Qtrait_loci{1,2}==idx - Ltot + N_of_Qtrait_loci) = [];

    state_of_Qtrait_loci{2,1}= state_of_Qtrait_loci{2,1}+1;

    state_of_Qtrait_loci{2,2} = [state_of_Qtrait_loci{2,2},idx- Ltot + N_of_Qtrait_loci];

end