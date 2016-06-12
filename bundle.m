%% Function bundles multiple simulations per file for when there are very many parameter sets (>>16)
function [fileNames] = bundle(parameter_matrix, rowCount, fileCount,dir_output)

partition = floor(rowCount/fileCount);

matPar = [];

for i=1:fileCount
    matPar = [matPar; partition];
end

parSum = sum(matPar);

if mod(rowCount, fileCount) > 0
    index = 1;
    while parSum < rowCount
        matPar(index) = matPar(index) + 1;
        index = index + 1;
        parSum = parSum + 1;
    end
end

fileNames = [];
rowStartIndex = 1;
cd(dir_output);

for i=1:fileCount
    rowEndIndex = rowStartIndex + matPar(i) - 1;
    subM = parameter_matrix(rowStartIndex:rowEndIndex, 1:14);    
    fileName = ['p' int2str(i)];    
    fileNames = strvcat(fileNames, fileName);
    dlmwrite(fileName, subM);
    rowStartIndex = rowEndIndex + 1;
end
