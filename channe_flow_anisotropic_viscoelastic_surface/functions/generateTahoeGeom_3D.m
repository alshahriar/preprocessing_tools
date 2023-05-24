function [xyzT,elmt,Nn,Ee,edgeList] = generateTahoeGeom_3D(x,y,z)
nxT = length(x);
nyT = length(y);
nzT = length(z);
Nn = nxT*nyT*nzT;
xyzT = zeros(Nn,3);

iN = 0;
for k = 1 : nzT
    for j = 1 : nyT
        for i = 1 : nxT
            iN = iN + 1;
            xyzT(iN,:) = [x(i) y(j) z(k)]; % coordinates of the tahoe grid points
        end
    end
end

% Connectivity
nXY = nxT*nyT;
Ee = (nzT-1)*(nyT-1)*(nxT-1);
elmt = zeros(Ee,8);
ie = 0;
for k = 1:(nzT-1)
    for j = 1:(nyT-1)
        for i = 1:(nxT-1)
            ie = ie + 1;
            elmt(ie,:) = [nXY*(k-1) + nxT*(j-1)+ i      nXY*(k-1) + nxT*(j-1)+i+1    nXY*(k-1) + nxT*j + i+1   nXY*(k-1) + nxT*j+ i ...
                nXY*(k-0) + nxT*(j-1)+ i      nXY*(k-0) + nxT*(j-1)+i+1    nXY*(k-0) + nxT*j + i+1   nXY*(k-0) + nxT*j+ i ];
        end
    end
end

edgeList = [];
% boundary nodes list - anticlock wise direction
% nPB = 2*nxT + 2*(nyT-2);
% edgeList = zeros(nPB,1);
% edgeList(1:nxT) = 1:nxT;
% edgeList(nxT+1:nxT + nyT-1) = (2:nyT)*nxT;
% edgeList(nxT+nyT-1 + 1 : nxT+nyT-1 + nxT-1) = nxT*nyT-1 : -1 : nxT*nyT-nxT+1;
% edgeList(nxT+nyT-1+nxT-1 + 1 : nPB) = (nyT-2:-1:1)*nxT + 1;

end