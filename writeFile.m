%%
function writeFile(outputFile,totalPoints,nElements)
fullFileName = fullfile('unstrucSurfBody1.dat');
fileID = fopen(fullFileName,'w');
fprintf(fileID,' TITLE = "Surface data"\n');
fprintf(fileID,'  variables = "x", "y", "z"\n');
fprintf(fileID,' zone N=    %d E=         %d\n',totalPoints,nElements);
fprintf(fileID,'  F=FEPOINT,ET=TRIANGLE\n');
fclose(fileID);
dlmwrite(fullFileName,outputFile(1:totalPoints,:),'-append', 'delimiter','\t', 'precision','%-.14f', 'newline', 'unix');
dlmwrite(fullFileName,outputFile(totalPoints+1:end,:),'-append', 'delimiter','\t', 'newline', 'unix');
end
