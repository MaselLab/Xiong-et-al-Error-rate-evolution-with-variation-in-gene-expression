function specificFileGen(dir_script,dir_ouput,dir_parameters,parNum)
%% load parameters
cd(dir_parameters);
fileID = fopen('specific_parameters_MC.txt');
fileMatrix = textscan(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 30);
fclose(fileID);

%% processing parameters matrix and packaging into smaller bundles
parameter_matrix = cell2mat(fileMatrix);
numRows = size(parameter_matrix, 1);
cd(dir_script);
fileNames = bundle(parameter_matrix, numRows, parNum,dir_ouput);
S = cellstr(fileNames);

%% run each bundle with a thread
cd(dir_script);
Pool=parpool('local',parNum);

[stream{1:parNum}]=RandStream.create('mrg32k3a','Seed','shuffle','Numstreams',parNum); % initilizing random number seed for each thread
spmd
    RandStream.setGlobalStream(stream{labindex});
end

numRows = size(fileNames, 1);
parfor i=1:numRows
    cd(dir_script);
    runBundle(S{i, 1},dir_script,dir_ouput,dir_parameters);
end

diary('off');
delete(Pool);