%% TahoeMarkers2FEMgridpointsBody1.dat
function writeFile2(relation)
%fullFileName = fullfile('C:\Users\alsha\Documents\result4', 'unstrucSurfBody1.dat');
fullFileName = fullfile('C:\Users\as17r\Dropbox\MatlabCodeForCreateGeomTecPlotFormat', 'TahoeMarkers2FEMgridpointsBody1.dat');
fileID = fopen(fullFileName,'w');
fclose(fileID);
dlmwrite(fullFileName,relation(1:end,1:4), 'delimiter','\t', 'precision','%-.14f', 'newline', 'unix');
end
