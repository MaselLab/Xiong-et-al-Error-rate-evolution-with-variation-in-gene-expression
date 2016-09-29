function MC_var_exp_plus_evo_biased(varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters and variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Loci Trans_exp_frq idx_of_ldel idx_of_lben Alpha Beta Grid_of_drho Grid_of_dalpha Grid_of_dbeta Grid_of_NewBeta PDF_of_drho PDF_of_dalpha PDF_of_dbeta PDF_of_NewBeta  %Psub_matrix_of_dalpha Psub_matrix_of_dbeta

%default initial conditions
Ltot = 100;
init_ldel =0;
POP_size = round(10^5.8);
log_of_rho = -2;
normal_distribution_sd = 2.25; %% this is changed manually. if you want, you can also modify formParameters to change it by input value
normal_distribution_mean = 11.5;
N_of_Qtrait = 1;
N_of_Qtrait_loci = 50;

%default mutation parameters
Mu_a = (23/9)*10^-9; %% per gene
Delta = log(10)*10^-2.5;
Sigma_f = 0.5;
Sigma_m = 0.5;
a = 750;
s = 20; 
Pdel = 0.4;
Pben = 0.1;
MutRate_of_loci = 3*10^-8;
MutRate_of_rho = 1*10^-6;
MutRate_of_alpha =3*10^-7;
MutRate_of_beta =3*10^-8;
new_optimal_trait =2;
SD_steady_state_beta = sqrt((Sigma_m/N_of_Qtrait_loci)^2/(1-((a-1)/a)^2));
bias_mut_rho=0;
SD_drho=0.2;

%default runtime parameters
outerLoopSteps = 105;
innerLoopSteps = 1000;
Range_drho = 2; % integral over -2 to 2
Range_of_dalpha = 5*Sigma_m/N_of_Qtrait_loci; % integrate over -5SD to 5SD
Range_of_dbeta = 5*Sigma_m/N_of_Qtrait_loci; 
Range_of_NewBeta = 5*SD_steady_state_beta;
GridNum = 2000;
GridNum_drho = GridNum;
burn_in = 100+1; % steps before changing optimal trait
recovery_flag1 = 0; %recovery of fitness
recovery_flag2 = 0; %recovery of Q trait
stop_flag =0;
recalc_mut_alphabeta = 1;

%modify the parameters if needed
[SD_drho,Mu_a, Delta, s, Pdel, Pben, MutRate_of_loci, MutRate_of_rho, MutRate_of_alpha, MutRate_of_beta, log_of_rho, outerLoopSteps, innerLoopSteps, Ltot, init_ldel ...
    , POP_size, normal_distribution_sd, normal_distribution_mean, N_of_Qtrait, N_of_Qtrait_loci, new_optimal_trait, a, bias_mut_rho,end_state, file_suffix, directory_of_source_code,directory_of_result] = formParameters(SD_drho,Mu_a, Delta, s, Pdel, Pben ...
    , MutRate_of_loci, MutRate_of_rho, MutRate_of_alpha, MutRate_of_beta, log_of_rho, outerLoopSteps, innerLoopSteps, Ltot, init_ldel, POP_size ...
    , normal_distribution_sd, normal_distribution_mean, N_of_Qtrait, N_of_Qtrait_loci, new_optimal_trait, a, bias_mut_rho,varargin, nargin);

% below are estimated sampling interval based on previous runs
% Here, I calculate the average simulation time (1/P_of_all) for 1000 mutation steps -- call it T1000.
% Then I divided T1000 by 10^7, and rounded the result.
% Rows correspond to ltot, 1000 to 200 from top to bottom.
% Columns correspond to popsize, 10^4 to 10^6 from left to right.
% Sampling with these intervals, different simulations should give similar numbers of data points.

interval_H=zeros(5,41);
interval_L=zeros(5,41);
entry1=linspace(1,41,21); % pop_size=10^2.0-10^6.0, increment=10^0.2
entry2=linspace(17,41,13); %pop_size=10^3.6-10^6.0, increment=10^0.2
entry3=[26,28,30,32];  %pop_size=10^(4.5,4.7,4.9,5.1)
switch normal_distribution_sd
    case 0
        interval_H(3,entry1)=[4.35,4.42,4.52,4.63,4.78,4.95,5.16,5.42,5.77,6.16,6.68,7.34,8.30,12.70,14.87,17.42,20.34,23.59,27.14,30.88,34.74];
        interval_H(3,entry3)=[11.76,13.73,16.09,18.84];  
        interval_L(3,entry1)=[4.35,4.42,4.52,4.63,4.78,4.95,5.16,5.42,5.77,6.16,6.68,7.34,8.14,9.04,9.98,10.91,20.34,23.60,27.16,30.88,34.72];
		interval_L(3,entry3)=[8.58,9.50,10.44,11.43];  
    case 0.6
        interval_H(3,entry2)=[5.76,6.16,6.67,7.33,8.14,12.57,14.89,17.41,20.27,23.54,27.15,30.78,34.71];  
        interval_L(3,entry1)=[5.77,6.17,6.70,7.33,8.15,9.03,9.97,10.97,20.31,23.58,27.11,30.90,34.68];	
    case 1.2
        interval_H(3,entry2)=[5.76,6.16,6.68,7.36,8.12,9.04,14.23,17.21,20.19,23.60,27.12,30.84,34.66];  
        interval_L(3,entry1)=[5.77,6.15,6.69,7.35,8.12,9.06,10.06,16.21,20.31,23.56,27.15,30.94,34.76];
    case 1.8
        interval_H(3,entry2)=[5.77,6.16,6.68,7.33,8.19,9.09,11.04,16.24,19.45,23.07,26.98,30.81,34.80];  
        interval_L(3,entry1)=[5.77,6.15,6.67,7.34,8.13,9.06,10.17,16.21,19.60,23.18,26.78,30.76,34.74];
    case 2.25
        interval_H(3,entry1)=[4.35,4.42,4.52,4.63,4.78,4.95,5.16,5.42,5.29,5.62,6.08,6.69,7.67,10.50,13.01,15.61,18.76,22.91,26.76,30.62,34.76];
        interval_H([1,2,4,5],entry2)=[  1.73,1.77,1.82,1.88,1.93,1.98,2.03,2.17,3.84,4.86,5.81,6.97,7.85;
                                        1.94,2.00,2.07,2.14,2.22,2.28,2.38,3.98,4.61,5.54,6.85,7.73,8.15;
                                        2.58,2.69,2.81,2.95,3.10,3.26,5.08,5.80,6.70,7.28,7.65,8.17,8.62;
                                        3.10,3.25,3.43,3.65,4.09,5.32,6.11,6.68,7.16,7.54,8.04,8.41,8.69];
        interval_H(3,entry3)=[8.63,9.71,13.32,16.70];  
        interval_L(3,entry1)=[4.35,4.42,4.52,4.63,4.78,4.95,5.16,5.42,5.77,6.16,6.69,7.35,8.16,9.13,10.50,14.89,18.27,22.02,26.13,30.26,34.44];
        interval_L([1,2,4,5],entry2)=[  1.73,1.78,1.83,1.88,1.93,1.98,2.03,2.20,3.81,4.95,6.03,6.74,7.85;
                                        1.93,2.00,2.07,2.14,2.21,2.28,2.38,3.51,4.80,5.85,6.78,7.47,8.21;
                                        2.59,2.69,2.81,2.94,3.07,3.35,4.43,5.28,6.53,7.30,7.86,8.29,8.62;
                                        3.10,3.25,3.43,3.66,4.14,5.57,5.95,6.59,7.13,7.65,7.98,8.34,8.68];
        interval_L(3,entry3)=[8.63,9.73,13.29,16.44];    
        if Pben/Pdel>=1
            interval_H(3,entry2)=[5.29,5.62,6.08,6.69,7.67,10.50,13.01,15.61,18.76,22.91,26.76,30.62,34.76];
            interval_L(3,entry2)=[5.28,5.62,6.07,6.70,7.71,10.42,12.86,15.37,19.28,22.72,26.83,30.28,34.82];     
        elseif Pben/Pdel>=0.01
            interval_H(3,entry2)=[7.13,7.63,8.50,9.71,11.10,12.93,15.11,17.81,21.18,24.75,28.97,34.67,41.40];
            interval_L(3,entry2)=[7.10,7.67,8.53,9.62,11.18,12.92,15.32,17.94,21.00,25.00,29.51,34.58,41.42]; 
        elseif Pben/Pdel>=0.001 % this condition was not tested in our study.
            interval_H(3,entry2)=[7.23,7.89,8.73,9.91,11.42,13.36,15.96,19.21,23.23,28.71,35.60,44.21,55.17];
            interval_L(3,entry2)=[7.28,7.86,8.66,9.84,11.37,13.25,15.95,19.37,23.41,28.97,35.84,44.38,55.04]; 
        end
        
        if bias_mut_rho~=0
            interval_H(3,entry2)=[5.77,6.16,6.68,7.35,8.17,9.14,10.56,14.88,18.36,22.07,26.06,30.21,34.37];
            interval_L(3,entry2)=[5.77,6.16,6.69,7.35,8.16,9.13,10.50,14.89,18.27,22.02,26.13,30.26,34.44];     
        end    
	
	if Delta< log(10)*10^-4.0
            interval_H(3,entry2)=[5.69,6.05,6.52,7.13,7.86,8.70,9.64,10.64,11.70,12.72,13.79,14.79,15.72];
            interval_L(3,entry2)=[5.69,6.05,6.52,7.13,7.86,8.70,9.65,10.64,11.68,12.73,13.77,14.81,15.85];
        elseif Delta< log(10)*10^-3.25
            interval_H(3,entry2)=[5.70,6.07,6.56,7.18,7.96,8.86,9.89,10.99,12.17,13.31,14.45,26.01,30.85];
            interval_L(3,entry2)=[5.70,6.07,6.56,7.18,7.95,8.86,9.89,11.01,12.14,13.30,14.49,17.78,30.63];
        elseif Delta< log(10)*10^-2.0
            interval_H(3,entry2)=[5.29,5.62,6.08,6.69,7.67,10.50,13.01,15.61,18.76,22.91,26.76,30.62,34.76];
            interval_L(3,entry2)=[5.28,5.62,6.07,6.70,7.71,10.42,12.86,15.37,19.28,22.72,26.83,30.28,34.82];
        elseif Delta< log(10)*10^-1.75
            interval_H(3,entry2)=[5.82,6.23,6.75,7.46,8.94,11.20,13.47,16.35,19.66,23.32,27.10,31.03,35.11];
            interval_L(3,entry2)=[5.82,6.22,6.74,7.47,8.86,11.18,13.59,16.23,19.70,23.36,27.15,31.10,35.14];
        elseif Delta< log(10)*10^-1.0
            interval_H(3,entry2)=[5.88,6.31,7.22,8.47,9.93,11.98,14.26,17.17,20.27,23.86,27.63,31.64,35.96];
            interval_L(3,entry2)=[5.88,6.32,7.24,8.40,9.94,11.87,14.20,17.20,20.41,23.89,27.65,31.68,35.98];
        elseif Delta< log(10)*10^-0.5
            interval_H(3,entry2)=[6.34,7.03,7.97,9.21,10.72,12.69,15.10,17.87,21.03,24.66,28.88,33.77,39.29];
            interval_L(3,entry2)=[6.34,7.03,7.98,9.20,10.80,12.76,15.15,17.87,21.05,24.65,28.85,33.79,39.42];
        end

    case 3.5
        interval_H(3,entry1)=[5.77,6.16,6.69,7.38,8.25,9.39,10.92,12.84,15.08,17.60,20.65,24.11,27.77]; 
        interval_H(3,entry3)=[8.78,10.09,11.87,13.74];
        interval_L(3,entry1)=[5.77,6.17,6.69,7.38,8.24,9.39,10.91,12.79,15.06,17.65,20.83,24.21,28.29];
        interval_L(3,entry3)=[8.74,10.10,11.86,13.90];  
end  

id_Ltot = 6-Ltot/200;

id_POPSIZE = round((log10(POP_size)-2.0)/0.1+1);

if init_ldel ==0; %high rho initial
    sample_interval = interval_H(id_Ltot,id_POPSIZE)*10^7;
else
    sample_interval = interval_L(id_Ltot,id_POPSIZE)*10^7;
end

time_to_sample = sample_interval;
count_N_sample = 1;
n_rec_step=1;
cd(directory_of_source_code);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(end_state,1)==1 % run simulation from start
    % loci, idx_of_lben, idx_of_ldel, 
    Loci = zeros(1,Ltot);
    idx = randperm(Ltot,init_ldel);
    Loci(idx) = 1; % deleterious cryptic sequence is 1
    idx_of_ldel =idx;
    idx_of_lben = find(Loci==0);
    sparse(Loci);
    N_of_ldel = init_ldel;
    state_of_Qtrait_loci = find_state_of_Qtrait_loci(Loci(Ltot-N_of_Qtrait_loci+1:Ltot));

    % Exp_frq and Trans_exp_frq
    Exp_frq = 2.^random('norm',normal_distribution_mean,normal_distribution_sd,Ltot,1); 
    Exp_frq = Exp_frq*(9000*Ltot/sum(Exp_frq));   
    Trans_exp_frq = Exp_frq';
    Tot_exp_frq = sum(Exp_frq);

    % Alpha, Beta, Oe
    Alpha = random('norm',0,SD_steady_state_beta,N_of_Qtrait_loci,N_of_Qtrait);
    Alpha = Alpha -mean(Alpha);
    Beta = random('norm',0,SD_steady_state_beta,N_of_Qtrait_loci,N_of_Qtrait);
    Beta = Beta - mean(Beta);
    optiaml_trait = zeros(1,N_of_Qtrait);  

    % rho
    rho = 10^log_of_rho;
    
else % use the end state of previous simulation to initialize everything
    % loci, idx_of_lben, idx_of_ldel,
    Loci = end_state(3,:);
    sparse(Loci);
    N_of_ldel = end_state(1,1);
    state_of_Qtrait_loci = find_state_of_Qtrait_loci(Loci(Ltot-N_of_Qtrait_loci+1:Ltot));
    idx_of_ldel = find(Loci==1);
    idx_of_lben = find(Loci==0);

    % % Exp_frq and Trans_exp_frq
    Trans_exp_frq=end_state(2,:);
    Exp_frq = Trans_exp_frq';
    Tot_exp_frq = sum(Exp_frq);

    % % Alpha, Beta, Oe
    Alpha = reshape(end_state(4,1:N_of_Qtrait*N_of_Qtrait_loci),N_of_Qtrait_loci,N_of_Qtrait);
    Beta = reshape(end_state(5,1:N_of_Qtrait*N_of_Qtrait_loci),N_of_Qtrait_loci,N_of_Qtrait);
    optiaml_trait = 2*ones(1,N_of_Qtrait);
    
    % rho
    rho = 10.^end_state(1,2);    
end

% drho grid and pdf
if SD_drho/0.2>1    
    GridNum_drho=GridNum_drho*SD_drho/0.2; % GridNum_drho is at least 2000. 0.2 is the default value of SD_drho
end
Range_drho=Range_drho*SD_drho/0.2;
Grid_of_drho = (linspace(-Range_drho,Range_drho,GridNum_drho))';
PD1 = ProbDistUnivParam('normal', [0 SD_drho]);
PDF_of_drho = pdf(PD1,Grid_of_drho);

% dalpha grid and pdf
Grid_of_dalpha = linspace(-Range_of_dalpha,Range_of_dalpha,GridNum);
PD2 = ProbDistUnivParam('normal', [0 Sigma_m/N_of_Qtrait_loci]);
PDF_of_dalpha = pdf(PD2,Grid_of_dalpha);

% dbeta grid and pdf
Grid_of_dbeta = linspace(-Range_of_dbeta,Range_of_dbeta,GridNum);
PD3 = ProbDistUnivParam('normal', [0 Sigma_m/N_of_Qtrait_loci]);
PDF_of_dbeta = pdf(PD3,Grid_of_dbeta);

% NewBeta grid and pdf
Grid_of_NewBeta = linspace(-Range_of_NewBeta,Range_of_NewBeta,GridNum);
PD4 = ProbDistUnivParam('normal', [0 SD_steady_state_beta]);
PDF_of_NewBeta = pdf(PD4,Grid_of_NewBeta);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% data storage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_waiting_time=0;

% record of alphas and betas
alphas_sampled_by_steps = zeros(outerLoopSteps+1,N_of_Qtrait*N_of_Qtrait_loci);
betas_sampled_by_steps = zeros(outerLoopSteps+1,N_of_Qtrait*N_of_Qtrait_loci);
alphas_sampled_by_time = zeros(outerLoopSteps*2,N_of_Qtrait*N_of_Qtrait_loci);
betas_sampled_by_time = zeros(outerLoopSteps*2,N_of_Qtrait*N_of_Qtrait_loci);
alphas_sampled_by_steps(1,:) = reshape(Alpha,1,N_of_Qtrait*N_of_Qtrait_loci);
betas_sampled_by_steps(1,:) = reshape(Beta,1,N_of_Qtrait*N_of_Qtrait_loci);
alphas_sampled_by_time(1,:) = reshape(Alpha,1,N_of_Qtrait*N_of_Qtrait_loci);
betas_sampled_by_time(1,:) = reshape(Beta,1,N_of_Qtrait*N_of_Qtrait_loci);

% record of ldel,rho and fitness 
record_sampled_by_steps=zeros(outerLoopSteps+1,6);
record_sampled_by_steps(1,3:4)=[init_ldel,log_of_rho];
record_sampled_by_time=zeros(outerLoopSteps*2,6);
record_sampled_by_time(1,3:4)=[init_ldel,log_of_rho];

% record the position of deleterious expression frequencies for each step
ldelIdx_sampled_by_steps = cell(outerLoopSteps+1, 2);
ldelIdx_sampled_by_steps{1,1} = 0;
ldelIdx_sampled_by_steps{1,2} = Loci;
ldelIdx_sampled_by_time = zeros(outerLoopSteps*2, Ltot);
ldelIdx_sampled_by_time(1,:) = Loci;


% matrix with all 6 fixation fluxes recorded per step
FixationFluxMatrix_sampled_by_steps = zeros(outerLoopSteps, 6);
FixationFluxMatrix_sampled_by_time = zeros(outerLoopSteps*2, 6);

% matrix that records 1/fixation_flux_of_all for each iteration of the outer loop
timeMatrix_sampled_by_steps = zeros(outerLoopSteps, 1);

% record cooption events
coopt_events = cell(outerLoopSteps*innerLoopSteps/10,6);
count_coopt_events = 1;

% record fitness and recovery time f
waiting_time_during_recovery = 0;
fitness_recovery_snapshot = zeros(7,1); % the elements are fitness before chaning Oe, fitness after changing Oe, fitness at 50% recovery, time in terms of 1/Pofall, and t in terms of steps, the step at which Oe changes, the step at which the criteria is met
Qtrait_recovery_snapshot = zeros(6,1); % the elements are Qtrait before chaning Oe, Qtrait at 50% recovery, time in terms of 1/Pofall, and t in terms of steps, the step at which Oe changes, the step at which the criteria is met

% output everything at every step during recovery
rec_alphas=zeros(500,N_of_Qtrait*N_of_Qtrait_loci);
rec_betas=zeros(500,N_of_Qtrait*N_of_Qtrait_loci);
rec_state=zeros(500,6);
rec_ldel_idx=zeros(500,Ltot);
rec_fixation_fluxes=zeros(500,6);

if size(end_state,1)~=1
    total_waiting_time=end_state(1,1);
    time_to_sample=total_waiting_time+sample_interval;
    record_sampled_by_time(1,2)=end_state(1,2);
    previous_steps=end_state(1,2);
    record_sampled_by_time(1,5)=end_state(1,5);
    record_sampled_by_time(1,6)=end_state(1,6);
    alphas_sampled_by_steps(1,:) = end_state(4,1:N_of_Qtrait*N_of_Qtrait_loci);
    betas_sampled_by_steps(1,:) = end_state(5,1:N_of_Qtrait*N_of_Qtrait_loci);
    alphas_sampled_by_time(1,:) = end_state(4,1:N_of_Qtrait*N_of_Qtrait_loci);
    betas_sampled_by_time(1,:) = end_state(5,1:N_of_Qtrait*N_of_Qtrait_loci);
else
    previous_steps=0;    
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:outerLoopSteps   % main iterations 
        if stop_flag ==1
            break;
        end    
       %% measure evolvability        
        if i == burn_in            
             [fitness_recovery_snapshot(1),Qtrait_recovery_snapshot(1),~,~] = calc_fitness(rho,s,Sigma_f,Delta,state_of_Qtrait_loci,N_of_Qtrait_loci,N_of_Qtrait,Tot_exp_frq,optiaml_trait,Exp_frq,Pben,Pdel);  % fitness and trait before optimal trait is changed           
             fitness_recovery_snapshot(6) = (i-1)*innerLoopSteps;  % when optimal trait is changed
             Qtrait_recovery_snapshot(5) = (i-1)*innerLoopSteps;
             waiting_time_during_recovery=0;             
             optiaml_trait = new_optimal_trait.*ones(1,N_of_Qtrait); 
	         recalc_mut_alphabeta = 1;           
             [fitness_recovery_snapshot(2),~,~,~] = calc_fitness(rho,s,Sigma_f,Delta,state_of_Qtrait_loci,N_of_Qtrait_loci,N_of_Qtrait,Tot_exp_frq,optiaml_trait,Exp_frq,Pben,Pdel);  % fitness after optiaml trait is changed            
        end
    
    for j = 1 : innerLoopSteps             
        if i >= burn_in && (recovery_flag1 == 0 || recovery_flag2 ==0)
            if recovery_flag1==0
                [fitness_recovery_snapshot(3),~,~,~] = calc_fitness(rho,s,Sigma_f,Delta,state_of_Qtrait_loci,N_of_Qtrait_loci,N_of_Qtrait,Tot_exp_frq,optiaml_trait,Exp_frq,Pben,Pdel) ;  % keep updating fitness until the loss of fitness is half recovered 
            end
            if recovery_flag2==0
                [~,Qtrait_recovery_snapshot(2),~,~]=calc_fitness(rho,s,Sigma_f,Delta,state_of_Qtrait_loci,N_of_Qtrait_loci,N_of_Qtrait,Tot_exp_frq,optiaml_trait,Exp_frq,Pben,Pdel) ; % keep updating trait until half optimal trait is reached
            end
            rec_alphas(n_rec_step,:)=reshape(Alpha,1,N_of_Qtrait*N_of_Qtrait_loci);
            rec_betas(n_rec_step,:)=reshape(Beta,1,N_of_Qtrait*N_of_Qtrait_loci);
            rec_ldel_idx(n_rec_step,:)=Loci;
            rec_state(n_rec_step,1)=waiting_time_during_recovery;
            rec_state(n_rec_step,2)=(i-1)*innerLoopSteps+j+previous_steps;
            rec_state(n_rec_step,3)=N_of_ldel;
            rec_state(n_rec_step,4)=log10(rho);
            [rec_state(n_rec_step,5),rec_state(n_rec_step,6),~,~]=calc_fitness(rho,s,Sigma_f,Delta,state_of_Qtrait_loci,N_of_Qtrait_loci,N_of_Qtrait,Tot_exp_frq,optiaml_trait,Exp_frq,Pben,Pdel);
            
            if  (2*fitness_recovery_snapshot(3)>=(fitness_recovery_snapshot(2)+ fitness_recovery_snapshot(1))) && recovery_flag1==0               
                fitness_recovery_snapshot(4) = waiting_time_during_recovery;
                fitness_recovery_snapshot(7) = (i-1)*innerLoopSteps+j;
                fitness_recovery_snapshot(5) = fitness_recovery_snapshot(7)-fitness_recovery_snapshot(6);
                recovery_flag1 = 1; %% the criteria has been met, so don't go through the above 
%                 stop_flag = 1; %% enable it if you want to stop the whole simulation at the end of the current inner loop after the fitness criteria is met
            end   
            
            if (2*Qtrait_recovery_snapshot(2)>=optiaml_trait)&&recovery_flag2==0
                recovery_flag2=1;
                Qtrait_recovery_snapshot(3)= waiting_time_during_recovery;
                Qtrait_recovery_snapshot(6)= (i-1)*innerLoopSteps+j;
                Qtrait_recovery_snapshot(4)=  Qtrait_recovery_snapshot(6)-Qtrait_recovery_snapshot(5);
            end
            
            if recovery_flag1==1 && recovery_flag2==1
                cd(directory_of_result); 
                dlmwrite(['fitness_recovery_snapshot' file_suffix],fitness_recovery_snapshot,'precision',8);
                dlmwrite(['Qtrait_recovery_snapshot' file_suffix],Qtrait_recovery_snapshot,'precision',8);
                cd(directory_of_source_code);
            end
        end        
          
       %% body of origin-fixation simulation
        Current_exp_of_ldel = Loci*Exp_frq;                
        [~,~,Current_Wm,trait_distance] = calc_fitness(rho,s,Sigma_f,Delta,state_of_Qtrait_loci,N_of_Qtrait_loci,N_of_Qtrait,Tot_exp_frq,optiaml_trait,Exp_frq,Pben,Pdel) ;
        
        % Calculate fixation flux of del2ben and ben2del
        [Pfix_of_ldel,Pfix_of_lben] = P_ben_and_del(POP_size,rho,s,Sigma_f,N_of_Qtrait_loci,Ltot,trait_distance,state_of_Qtrait_loci,Current_Wm,Current_exp_of_ldel,N_of_ldel,Tot_exp_frq,Pben,Pdel);        
        P_of_ldel = POP_size*MutRate_of_loci*Pdel*sum(Pfix_of_ldel);         
        P_of_lben = POP_size*MutRate_of_loci*Pben*sum(Pfix_of_lben); 

        % Calculate fixation flux of rho mutation
        [P_of_drho, max_Pfix_of_drho] = P_drho(Delta,MutRate_of_rho,POP_size,rho,s,Sigma_f,N_of_Qtrait_loci,trait_distance,state_of_Qtrait_loci,Current_exp_of_ldel,Tot_exp_frq,Range_drho,GridNum_drho,Pben,Pdel); %% replaced by new function P_drho
        P_of_drho =POP_size*P_of_drho;
	
        % Calculate fixation flux of alpha and beta mutation    
        if recalc_mut_alphabeta ~=0  %% if mutation of cryptic state is outside quantitative loci, don't recalculate
            [P_of_dalpha,P_of_dbeta, max_Pfix_of_dalpha, max_Pfix_of_dbeta] = P_alpha_beta_biased(MutRate_of_alpha,MutRate_of_beta,POP_size,rho,Sigma_f,...
                                                                              N_of_Qtrait_loci,trait_distance,state_of_Qtrait_loci,N_of_Qtrait,Range_of_dalpha,Range_of_dbeta,GridNum,a); %% new function plugged in
      	    P_of_dalpha = POP_size*P_of_dalpha;
            P_of_dbeta = POP_size*P_of_dbeta;
        end
                  
        % calculate fixation flux of cooption
        [fixation_flux_of_cooption, fixation_flux_ben, fixation_flux_ben_per_locus, max_Pfix_of_NewBeta_ben, fixation_flux_del, fixation_flux_del_per_locus] = ...
        P_coopt(rho,Sigma_f,state_of_Qtrait_loci,trait_distance,POP_size, Pdel,Pben,Current_exp_of_ldel,Current_Wm,Range_of_NewBeta,Ltot, N_of_Qtrait_loci,GridNum,N_of_Qtrait,s,Tot_exp_frq); %% new function, it returns a matrix of the probability of cooption happening at every locus in the K_set;         
        P_of_coopt = POP_size*Mu_a*fixation_flux_of_cooption;  
        
        
        % pick up mutations at each step
        P_of_all = P_of_ldel+ P_of_lben + P_of_drho + P_of_dalpha + P_of_dbeta + P_of_coopt;
        Mut = datasample([1,2,3,4,5,6],1,'Weight',[P_of_ldel/P_of_all,P_of_drho/P_of_all,P_of_lben/P_of_all,P_of_dalpha/P_of_all,P_of_dbeta/P_of_all,P_of_coopt/P_of_all]); % 1 for lben to ldel, 2 for mutating rho, 3 for ldel to lben

        % generate mutations
        recalc_mut_alphabeta = 1;
        
        switch Mut
            case 2 % rho
                rho = muta_rho(Delta,POP_size,rho,s,Sigma_f,Current_exp_of_ldel,Tot_exp_frq,state_of_Qtrait_loci,trait_distance,N_of_Qtrait_loci,max_Pfix_of_drho,Pben,Pdel,bias_mut_rho,SD_drho); %% sampling and mutaion are now merged in one function. rewrote the code so that muta_rho, muta_alpha and muta_beta each contains a rejection method module.                
            case 4 % alpha
                muta_alpha_beta_biased(rho,POP_size,Sigma_f,Sigma_m,trait_distance,N_of_Qtrait_loci,N_of_Qtrait,max_Pfix_of_dalpha,a,Ltot,1); %% same as above                               
            case 5 % beta
                muta_alpha_beta_biased(rho,POP_size,Sigma_f,Sigma_m,trait_distance,N_of_Qtrait_loci,N_of_Qtrait,max_Pfix_of_dbeta,a,Ltot,2); %% same as above            	                
            case 6 % cooption
                [idx_of_coopted_locus,N_of_ldel,state_of_Qtrait_loci] = muta_coopt(state_of_Qtrait_loci, SD_steady_state_beta,rho,trait_distance, N_of_ldel, fixation_flux_ben,...
                                                                 fixation_flux_ben_per_locus, max_Pfix_of_NewBeta_ben, fixation_flux_del, fixation_flux_del_per_locus,Sigma_f,POP_size,Ltot,N_of_Qtrait_loci); 
                coopt_events{count_coopt_events,1} = (i-1)*innerLoopSteps+j;
                coopt_events{count_coopt_events,2} = rho;
                coopt_events{count_coopt_events,3} = Loci(Ltot-N_of_Qtrait_loci+1:Ltot);
                coopt_events{count_coopt_events,4} = idx_of_coopted_locus + Ltot-N_of_Qtrait_loci;
                coopt_events{count_coopt_events,5} = Alpha(idx_of_coopted_locus,:);
                coopt_events{count_coopt_events,6} = Beta(idx_of_coopted_locus,:);
                count_coopt_events = count_coopt_events+1;                
            case 3 % ldel to lben
                [N_of_ldel,state_of_Qtrait_loci,recalc_mut_alphabeta] = mutate_ldel_to_lben(Ltot,N_of_Qtrait_loci, Pfix_of_lben, N_of_ldel,state_of_Qtrait_loci);                         
            case 1 % lben to ldel
                [N_of_ldel,state_of_Qtrait_loci,recalc_mut_alphabeta] = mutate_lben_to_ldel(Ltot,N_of_Qtrait_loci, Pfix_of_ldel, N_of_ldel,state_of_Qtrait_loci);               
        end
       
       %% extract data
        total_waiting_time=total_waiting_time+round(1/P_of_all);
        waiting_time_during_recovery = waiting_time_during_recovery+round(1/P_of_all);

        if total_waiting_time >= time_to_sample                       
            alphas_sampled_by_time(count_N_sample+1,:) = reshape(Alpha,1,N_of_Qtrait*N_of_Qtrait_loci);
            betas_sampled_by_time(count_N_sample+1,:) = reshape(Beta,1,N_of_Qtrait*N_of_Qtrait_loci);            
            ldelIdx_sampled_by_time(count_N_sample+1,:) = Loci;
            record_sampled_by_time(count_N_sample+1,1) = total_waiting_time;
            record_sampled_by_time(count_N_sample+1,2) = (i-1)*innerLoopSteps+j+previous_steps;
            record_sampled_by_time(count_N_sample+1,3) = N_of_ldel;
            record_sampled_by_time(count_N_sample+1,4) = log10(rho);
            [record_sampled_by_time(count_N_sample+1,5),record_sampled_by_time(count_N_sample+1,6),~,~] = calc_fitness(rho,s,Sigma_f,Delta,state_of_Qtrait_loci,N_of_Qtrait_loci,N_of_Qtrait,Tot_exp_frq,optiaml_trait,Exp_frq,Pben,Pdel);
            FixationFluxMatrix_sampled_by_time(count_N_sample, 1) = P_of_ldel;
            FixationFluxMatrix_sampled_by_time(count_N_sample, 2) = P_of_lben;
            FixationFluxMatrix_sampled_by_time(count_N_sample, 3) = P_of_drho;
            FixationFluxMatrix_sampled_by_time(count_N_sample, 4) = P_of_dalpha;
            FixationFluxMatrix_sampled_by_time(count_N_sample, 5) = P_of_dbeta;
            FixationFluxMatrix_sampled_by_time(count_N_sample, 6) = P_of_coopt;            
            time_to_sample = time_to_sample+ sample_interval;
            count_N_sample = count_N_sample +1;
        end
        
        if i >= burn_in && (recovery_flag1 == 0 || recovery_flag2==0)      % record fixation fluxes during recovery
            rec_fixation_fluxes(n_rec_step,1)= P_of_ldel;
            rec_fixation_fluxes(n_rec_step,2)= P_of_lben;
            rec_fixation_fluxes(n_rec_step,3)= P_of_drho;
            rec_fixation_fluxes(n_rec_step,4)= P_of_dalpha;
            rec_fixation_fluxes(n_rec_step,5)= P_of_dbeta;
            rec_fixation_fluxes(n_rec_step,6)= P_of_coopt;
            n_rec_step=n_rec_step+1;
        end       
      
        if j==innerLoopSteps-1                                             % record position of loci, number of ldel, and rho for the
            ldelIdx_sampled_by_steps{i+1, 1} = i*innerLoopSteps;           % 100th state (which appears whenever the inner loop is at the
            ldelIdx_sampled_by_steps{i+1, 2} = Loci;                       % 99th iteration)
            record_sampled_by_steps(i+1,1)=total_waiting_time;
            record_sampled_by_steps(i+1,2) = i;
            record_sampled_by_steps(i+1,3)= N_of_ldel;
            record_sampled_by_steps(i+1,4)= log10(rho);
            [record_sampled_by_steps(i+1,5),record_sampled_by_steps(i+1,6),~,~]= calc_fitness(rho,s,Sigma_f,Delta,state_of_Qtrait_loci,N_of_Qtrait_loci,N_of_Qtrait,Tot_exp_frq,optiaml_trait,Exp_frq,Pben,Pdel);
            alphas_sampled_by_steps(i+1,:) = reshape(Alpha,1,N_of_Qtrait*N_of_Qtrait_loci);
            betas_sampled_by_steps(i+1,:) = reshape(Beta,1,N_of_Qtrait*N_of_Qtrait_loci);   
        end

    end % end of inner loop
    
    % record fixation fluxes
    FixationFluxMatrix_sampled_by_steps(i, 1) = P_of_ldel;
    FixationFluxMatrix_sampled_by_steps(i, 2) = P_of_lben;
    FixationFluxMatrix_sampled_by_steps(i, 3) = P_of_drho;
    FixationFluxMatrix_sampled_by_steps(i, 4) = P_of_dalpha;
    FixationFluxMatrix_sampled_by_steps(i, 5) = P_of_dbeta;
    FixationFluxMatrix_sampled_by_steps(i, 6) = P_of_coopt;
    
    timeMatrix_sampled_by_steps(i) = 1/P_of_all;
    
end % end of outerloop

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% data packaging and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphas_sampled_by_time(count_N_sample+1:end,:) = [];
betas_sampled_by_time(count_N_sample+1:end,:) = [];            
ldelIdx_sampled_by_time(count_N_sample+1:end,:) = [];
record_sampled_by_time(count_N_sample+1:end,:)=[];
FixationFluxMatrix_sampled_by_time(count_N_sample:end,:) = [];
rec_alphas(n_rec_step+1:end,:)=[]; 
rec_betas(n_rec_step+1:end,:)=[];
rec_state(n_rec_step+1:end,:)=[];
rec_ldel_idx(n_rec_step+1:end,:)=[];
rec_fixation_fluxes(n_rec_step+1:end,:)=[];

cd(directory_of_result);
dlmwrite(['rec_alphas' file_suffix],rec_alphas,'precision',8);
dlmwrite(['rec_betas' file_suffix],rec_betas,'precision',8);
dlmwrite(['rec_state' file_suffix],rec_state,'precision',8);
dlmwrite(['rec_probabilities' file_suffix],rec_fixation_fluxes,'precision',8);
dlmwrite(['rec_ldel_idx' file_suffix],rec_ldel_idx,'precision',8);
dlmwrite(['final_state' file_suffix],record_sampled_by_steps,'precision',8);
dlmwrite(['final_state_sampled_by_time' file_suffix],record_sampled_by_time,'precision',8);
dlmwrite(['Exp_frq' file_suffix],Exp_frq,'precision',8);
dlmwrite(['Alphas' file_suffix], alphas_sampled_by_steps, 'precision',8);
dlmwrite(['Alphas_sampled_by_time' file_suffix], alphas_sampled_by_time, 'precision',8);
dlmwrite(['Betas' file_suffix], betas_sampled_by_steps, 'precision',8);
dlmwrite(['Betas_sampled_by_time' file_suffix], betas_sampled_by_time, 'precision',8);
dlmwrite(['1divP_of_all' file_suffix], timeMatrix_sampled_by_steps,'precision',8);
dlmwrite(['probabilities' file_suffix], FixationFluxMatrix_sampled_by_steps,'precision',8);
dlmwrite(['probabilities_sampled_by_time' file_suffix], FixationFluxMatrix_sampled_by_time,'precision',8);
dlmwrite(['coopt_events' file_suffix], coopt_events,'precision',8);
ldelIdx_sampled_by_steps = cell2mat(ldelIdx_sampled_by_steps);
dlmwrite(['ldel_idx_tot' file_suffix], ldelIdx_sampled_by_steps,'precision',8);
dlmwrite(['ldel_idx_tot_sampled_by_time' file_suffix], ldelIdx_sampled_by_time,'precision',8);
cd(directory_of_source_code);

clearvars -global;

