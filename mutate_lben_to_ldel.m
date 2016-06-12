function[N_of_ldel,state_of_Qtrait_loci,recalc_Wq] = mutate_lben_to_ldel(Ltot,N_of_Qtrait_loci, Pfix_of_ldel, N_of_ldel, state_of_Qtrait_loci)
global Loci idx_of_ldel idx_of_lben

Pfix_of_ldel=Pfix_of_ldel/sum(Pfix_of_ldel);

idx = datasample(idx_of_lben,1,'Weight',Pfix_of_ldel);

recalc_Wq = 1;

if idx < Ltot-N_of_Qtrait_loci+1   %% if the mutation is outside the K set, no need to recalculate Wtrait in the next round 
    
    recalc_Wq = 0;
    
end

Loci(idx) = 1;

idx_of_lben(idx_of_lben==idx)=[];

idx_of_ldel = [idx_of_ldel,idx];

N_of_ldel = N_of_ldel+1;

if idx > Ltot - N_of_Qtrait_loci

    state_of_Qtrait_loci{2,1}= state_of_Qtrait_loci{2,1}-1;

    state_of_Qtrait_loci{2,2}(state_of_Qtrait_loci{2,2}==idx - Ltot + N_of_Qtrait_loci) = [];

    state_of_Qtrait_loci{1,1}= state_of_Qtrait_loci{1,1}+1;

    state_of_Qtrait_loci{1,2} = [state_of_Qtrait_loci{1,2},idx- Ltot + N_of_Qtrait_loci];

end