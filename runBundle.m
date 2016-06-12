function runBundle(simu_conditions, dir_script,dir_output,dir_parameters)

cd(dir_output);
simu_con_matrix=csvread(simu_conditions);

N_simu_conditions = size(simu_con_matrix, 1);
cd(dir_script);

for i=1:N_simu_conditions
    MC_var_exp_plus_evo_biased(simu_conditions, i,dir_script,dir_output,dir_parameters); 
end
