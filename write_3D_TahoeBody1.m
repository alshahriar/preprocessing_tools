function write_3D_TahoeBody1(xyz,elmt,Nn,Ee)
%TITLE = "Example: FE-Volume Brick Data"
%VARIABLES = "X", "Y", "Z", "Temperature"
%ZONE NODES=27, ELEMENTS=6, DATAPACKING=POINT, ZONETYPE=FEBRICK
% in Tahoe format - grid file for Tahoe
fid = fopen('TahoeBody1.dat','w');
fprintf(fid,'TITLE = "TAHOE FEM Grid data Brick"\n');
fprintf(fid,'VARIABLES = "X", "Y", "Z" \n');
fprintf(fid,'ZONE NODES=%d, ELEMENTS=%d\n',Nn,Ee);
fprintf(fid,',Datapacking=point,Zonetype=FEBRICK\n');

for iN = 1:Nn
    fprintf(fid,'%15e %15e %15e\n',xyz(iN,1),xyz(iN,2),xyz(iN,3));
end

for ie = 1:Ee
    fprintf(fid,'%12d %12d %12d %12d %12d %12d %12d %12d\n',elmt(ie,1),elmt(ie,2),elmt(ie,3),elmt(ie,4),elmt(ie,5),elmt(ie,6),elmt(ie,7),elmt(ie,8));
end

fclose(fid);
end
