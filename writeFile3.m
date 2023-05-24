function writeFile3(relation)
%fullFileName = fullfile('C:\Users\alsha\Documents\result4', 'unstrucSurfBody1.dat');
fullFileName = fullfile('C:\Users\as17r\Dropbox\MatlabCodeForCreateGeomTecPlotFormat', 'LeftRight.dat');
fileID = fopen(fullFileName,'w');
fclose(fileID);
dlmwrite(fullFileName,relation(1:end,1:3), 'delimiter','\t', 'newline', 'unix');
end
