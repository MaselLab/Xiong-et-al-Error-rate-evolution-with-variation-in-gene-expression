function[Miu, Delta, s, Pdel, Pben, MutRate_of_loci, MutRate_of_rho, MutRate_of_alpha, MutRate_of_beta, log_of_rho, outerLoopSteps, innerLoopSteps, Ltot, init_ldel, POP_size, normal_distribution_sd, normal_distribution_mean, N_of_Qtrait, N_of_Qtrait_loci, delta_oe, a, bias_mut_rho, end_state,file_suffix,dir_script,dir_output] ...
    = formParameters(Miu, Delta, s, Pdel, Pben, MutRate_of_loci, MutRate_of_rho, MutRate_of_alpha, MutRate_of_beta, log_of_rho, outerLoopSteps, innerLoopSteps, Ltot, init_ldel, POP_size, normal_distribution_sd, normal_distribution_mean, N_of_Qtrait, N_of_Qtrait_loci, delta_oe, a, bias_mut_rho, arguments, numvars)
%% if simulation condition is provided, take it
use_input = false; % default value
if numvars ==5
    use_input = true;
    simu_conditions = arguments{1};   
    dir_script = arguments{3};
    dir_output = arguments{4};
    dir_parameters =arguments{5};
else
    dir_script = arguments{1};
    dir_output = arguments{2};
    dir_parameters =arguments{3};
end    

if use_input == true    
    cd(dir_output);
    matrixOfParameters = csvread(simu_conditions);   
    rowToRun = arguments{2};
    
    param1 = matrixOfParameters(rowToRun, 1);
    val1 = matrixOfParameters(rowToRun, 2);
    param2 = matrixOfParameters(rowToRun, 3);
    val2 = matrixOfParameters(rowToRun, 4);
    param3 = matrixOfParameters(rowToRun, 5);
    val3 = matrixOfParameters(rowToRun, 6);
    param4 = matrixOfParameters(rowToRun, 7);
    val4 = matrixOfParameters(rowToRun, 8);
    param5 = matrixOfParameters(rowToRun, 9);
    val5 = matrixOfParameters(rowToRun, 10);
    param6 = matrixOfParameters(rowToRun, 11);
    val6 = matrixOfParameters(rowToRun, 12);
    param7 = matrixOfParameters(rowToRun, 13);
    val7 = matrixOfParameters(rowToRun, 14);

    matParam = [param1, param2, param3, param4, param5, param6, param7];
    matVal = [val1, val2, val3, val4, val5, val6, val7];
    allowedParameterAmount = 7; 

    rep_id=-1;	%By default, there is one replicate for each condition, therefore do not distinguish replicates in filenames. 
    cont=0; %By default, simulations always start from the beginning, rather than continue from a previous simulation. 
        
    for i=1:allowedParameterAmount % switches the specified parameters with their new values
        if matParam(i) > 0
            paramTemp = matParam(i);
            valTemp = matVal(i);
            switch paramTemp                
                case 1
                    Miu = valTemp;               
                case 2
                    Delta = valTemp;
                case 3
                    s = valTemp;
                case 4
                    Pdel = valTemp;
                case 5
                    Pben = valTemp;
                case 6
                    MutRate_of_loci = valTemp;
                case 7
                    MutRate_of_rho = valTemp;
                case 8
                    init_condition = valTemp;
                case 9
                    N_of_Qtrait_loci = valTemp;
                case 10
                    bias_mut_rho = valTemp;
                case 11
                    Ltot = valTemp;                
                case 12
                    POP_size = round(10^valTemp);
                case 13
                    MutRate_of_alpha = valTemp;
                case 14
                    MutRate_of_beta = valTemp;
                case 15
                    rep_id = valTemp;  
                case 16
                    a=valTemp;
                case 17
                    cont=1;
            end
        end
    end
end 

if init_condition==1 % high-rho initial condition
	init_ldel=0;
	log_of_rho=-1.0;
else 
	init_ldel=round(Ltot*Pdel/(Pben+Pdel));
	log_of_rho=-5.0;
end

%% Make a report of important parameters
strMua = num2str(Miu);
strSigma = num2str(Delta);
strS = num2str(s);
strPdel = num2str(Pdel);
strPben = num2str(Pben);
strMutRateLoci = num2str(MutRate_of_loci);
strMutRateRho = num2str(MutRate_of_rho);
strSteps = int2str(outerLoopSteps);
strSkip = int2str(innerLoopSteps);
strLdel = int2str(init_ldel);
strPop = int2str(POP_size);
strLtot = int2str(Ltot);
strNormDistSD = num2str(normal_distribution_sd);
strNormDistMean = int2str(normal_distribution_mean);
strRho = num2str(log_of_rho);
strMutRateAlpha = num2str(MutRate_of_alpha);
strMutRateBeta = num2str(MutRate_of_beta);
strQtrait = num2str(N_of_Qtrait);
strQtraitLoci = num2str(N_of_Qtrait_loci);
strDeltaOE = num2str(delta_oe);
strA = num2str(a);
strbias_rho = num2str(bias_mut_rho);
if rep_id ~=-1
        rep_id = num2str(rep_id);
	strRepID = strcat('_rep',rep_id);
else 
	strRepID = '';
end
vectorOfParameters = {  ['Mu_a: ' strMua];
                        ['Delta: ' strSigma];
                        ['s: ' strS];
                        ['Pdel: ' strPdel];
                        ['Pben: ' strPben];
                        ['MutRate_of_loci: ' strMutRateLoci];
                        ['MutRate_of_rho: ' strMutRateRho];
                        ['MutRate_of_alpha: ' strMutRateAlpha];
                        ['MutRate_of_beta: ' strMutRateBeta];
                        ['Outer loop itr: ' strSteps];
                        ['Inner loop itr: ' strSkip];
                        ['initial ldel: ' strLdel];
                        ['N: ' strPop];
                        ['Ltot: ' strLtot];
                        ['Distribution SD: ' strNormDistSD];
                        ['Distribution Mean: ' strNormDistMean];
                        ['Log10(rho): ' strRho];
                        ['N_of_Qtrait: ' strQtrait];
                        ['N_of_Qtrait_loci: ' strQtraitLoci];
                        ['delta_oe: ' strDeltaOE];
                        ['a: ' strA];
                        ['mut_bias_rho: ', strbias_rho]};
                    
% make strings from a few of the parameters for the filename suffix
file_suffix = strcat(strPop, 'x', strLtot, strRepID);

cd(dir_output);
fileID = fopen(['parameters' file_suffix], 'w');
formatSpec = '%s\n';
for i=1:size(vectorOfParameters, 1);
    fprintf(fileID, formatSpec, vectorOfParameters{i});
end
fclose(fileID);

%% if continue simulation using end_state from previous simulations
if cont==1      
    cd(dir_parameters);
    filename=strcat('End_state_',strPop,'x',strRepID);
    f_Estate=dir(filename);
    end_state=load(f_Estate.name);
else
    end_state=0;
end
