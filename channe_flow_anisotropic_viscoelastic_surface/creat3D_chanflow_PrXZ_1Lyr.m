% this is for Tahoe - Double layer of triangles

clear all
geometry = 'rectangle';

if(strcmp(geometry,'rectangle'))
    
    if(0) % Debug setup
        delta = 1;
        xOffset = 1;
        zOffset = 1;
        xs = -xOffset; ys = 0; zs = -zOffset; % left bottom corner
        lengthX = 2;
        lengthY = 2;
        lengthZ = 2;
        %lengthY = 1.001; % outBItriangle.err
        %lengthY = 1; % outBItriangle.err
        xe = xs + lengthX + 2*xOffset;
        ye = ys + lengthY;
        nx = 5;
        ny = 3;
        leftSide = 0; rightSide = 2;
        frontSide = 0; backSide = 2;
    elseif(0)  % Debug setup
        delta = 1;
        xOffset = 1;
        zOffset = 1;
        xs = -xOffset; ys = 0; zs = -zOffset; % left bottom corner
        lengthX = 3;
        lengthY = 2;
        lengthZ = 2;
        %lengthY = 1.001; % outBItriangle.err
        %lengthY = 1; % outBItriangle.err
        xe = xs + lengthX + 2*xOffset;
        ye = ys + lengthY;
        nx = 6;
        ny = 3;
        leftSide = 0; rightSide = 3;
        frontSide = 0; backSide = 2;
    elseif(0)  % Debug setup
        delta = 1;
        xOffset = 1;
        zOffset = 1;
        xs = -xOffset; ys = 0; zs = -zOffset; % left bottom corner
        lengthX = 4;
        lengthY = 2;
        lengthZ = 2;
        %lengthY = 1.001; % outBItriangle.err
        %lengthY = 1; % outBItriangle.err
        xe = xs + lengthX + 2*xOffset;
        ye = ys + lengthY;
        nx = 7;
        ny = 3;
        leftSide = 0; rightSide = 4;
        frontSide = 0; backSide = 2;
    end
    
    if(1) % Setup-1: Channel flow
        delta = 1;
        xOffset = 0.25;
        zOffset = xOffset;
        xs = -xOffset; ys = 0; zs = -zOffset; % left bottom corner
        lengthX = 2*3.1416;
        lengthY = 2;
        lengthZ = 3.1416;
        %lengthY = 1.001; % outBItriangle.err
        %lengthY = 1; % outBItriangle.err
        xe = xs + lengthX + 2*xOffset;
        ye = ys + lengthY;
        nx = 100;
        ny = 15;
        leftSide = 0; rightSide = 2*3.1416;
        frontSide = 0; backSide = 3.1416;
    elseif(0) % Setup-2: not working
        delta = 1;
        xOffset = 0.4;
        zOffset = xOffset;
        xs = -xOffset; ys = 0; zs = -zOffset; % left bottom corner
        lengthX = 2*pi;
        lengthY = 2;
        lengthZ = pi;
        %lengthY = 1.001; % outBItriangle.err
        %lengthY = 1; % outBItriangle.err
        xe = xs + lengthX + 2*xOffset;
        ye = ys + lengthY;
        nx = 20;
        ny = 3;
        leftSide = 0; rightSide = 2*pi;
        frontSide = 0; backSide = pi;
    end
    [newX, newY, xOffset, x, y] = get_nexX_newY(xs,xe,nx,ys,ye,ny,xOffset,leftSide,rightSide);
end
dx = x(2) - x(1);
disp(['Boundary layer thickness = ',num2str(delta)]);
%% 
%       pointsP1(:,1) ----> x coordinate
%       pointsP1(:,2) ----> y coordinate

nP_newX = length(newX);
dz = dx; % same as dx
[newTempZ,z] = get_newTempZ(zs,lengthZ,zOffset,frontSide,backSide,dz);
nz = length(z);
nPlaneInZ = length(newTempZ);
for i = 1:nPlaneInZ
    pointsP1((i-1)*nP_newX+1:i*nP_newX,:)=[newX' newY' newX*0+newTempZ(i)];
end
newZ = pointsP1(:,3)';

%% generating the relation of the nodes for the elements and writing the file unstrucSurfBody1.dat
nPoints = length(newX);
nElements = nPoints*2*(nPlaneInZ-1);

for kL = 1:nPlaneInZ-1
    coeff = (kL-1)*nPoints;
    for i=1:nPoints-1
        elem(2*coeff + 2*i-1,:)=[coeff+i,coeff+i+1, coeff+nPoints+i]; %setting index for each elements
        elem(2*coeff + 2*i,:)=[coeff+i+1,  coeff+nPoints+i+1, coeff+nPoints+i];
    end
    %For the last element
    elem(2*coeff+2*nPoints-1,:)=[coeff+nPoints,coeff+1, coeff+nPoints+nPoints];
    elem(2*coeff+2*nPoints,:)=[coeff+1, coeff+nPoints+1,coeff+nPoints+nPoints];
end
%
%
% %1st Layer of the Elements
% for i=1:nPoints-1
%     elem(2*i-1,:)=[i,i+1, nPoints+i]; %setting index for each elements
%     elem(2*i,:)=[i+1,  nPoints+i+1, nPoints+i];
% end
% %For the last element
% elem(2*nPoints-1,:)=[nPoints,1, nPoints+nPoints];
% elem(2*nPoints,:)=[1, nPoints+1,nPoints+nPoints];
%
% %2nd Layer of the Elements
% for i=1:nPoints-1
%     elem(2*nPoints+2*i-1,:)=[nPoints+i,nPoints+i+1, nPoints+nPoints+i]; %setting index for each elements
%     elem(2*nPoints+2*i,:)=[nPoints+i+1,  nPoints+nPoints+i+1, nPoints+nPoints+i];
% end
% %For the last element
% elem(2*nPoints+2*nPoints-1,:)=[nPoints+nPoints,nPoints+1, nPoints+nPoints+nPoints];
% elem(2*nPoints+2*nPoints,:)=[nPoints+1, nPoints+nPoints+1,nPoints+nPoints+nPoints];
%
outputFile=[pointsP1;elem];
writeFile(outputFile,nPlaneInZ*nPoints,nElements); % writing the file unstrucSurfBody1.dat with double layer

%% Making Grids for Tahoe solver
% some part of the body will be flexible and some are rigid
yTop=max(y);
yBottom=min(y);

for i=1:nx
    if(x(i)>leftSide && x(i)<(leftSide+dx))
        leftMost=i;
    elseif (abs(x(i)-leftSide)<1e-8)
        leftMost=i;break;
    end
end

for i=1:nx
    if(x(i)<rightSide && x(i)>(rightSide-dx))
        rightMost=i;% index of another point
    elseif (abs(x(i)-rightSide)<1e-8)
        rightMost=i;break;
    end
end

for k=1:nz
    if(z(k)>frontSide && z(k)<(frontSide+dz))
        frontMost=k; %index of P2 point
    elseif (abs(z(k)-frontSide)<1e-8)
        frontMost=k;break;
    end
end

for k=1:nz-1
    if(z(k)<backSide && z(k+1)>(backSide))
        backMost=k;% index of another point
    elseif(abs(z(k)-backSide)<1e-8)
        backMost=k;break;
    end
end

xTahoe = x(leftMost:rightMost);
yTahoe = y;
zTahoe = z(frontMost:backMost);

nxT = length(xTahoe);
nyT = length(yTahoe);
nzT = length(zTahoe);

Nn = nxT*nyT*nzT;

pointTemp = zeros(Nn,3);

iN = 0;
for k = 1 : nzT
    for j = 1 : nyT
        for i = 1 : nxT
            iN = iN + 1;
            pointTemp(iN,:) = [xTahoe(i) yTahoe(j) zTahoe(k)]; % coordinates of the tahoe grid points
        end
    end
end

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
            %            elmt(ie,:) = [nxT*(j-1)+i     nxT*(j-1)+i+1    nxT*j+i+1   nxT*j+i];
        end
    end
end


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
    fprintf(fid,'%15e %15e %15e\n',pointTemp(iN,1),pointTemp(iN,2),pointTemp(iN,3));
end

for ie = 1:Ee
    fprintf(fid,'%12d %12d %12d %12d %12d %12d %12d %12d\n',elmt(ie,1),elmt(ie,2),elmt(ie,3),elmt(ie,4),elmt(ie,5),elmt(ie,6),elmt(ie,7),elmt(ie,8));
end

fclose(fid);



%% Connecting points between IBM body grid and Tahoe grid

nPointsForUnstructured = length(pointsP1(:,1));

fid = fopen('TahoeMarkers2FEMgridpointsBody1.dat','w');

for iN = 1:nPointsForUnstructured
    for k = 1:nzT
        for j = 1: nyT
            for i = 1: nxT
                if(abs(pointsP1(iN,1)-xTahoe(i))<1e-8 && abs(pointsP1(iN,2)-yTahoe(j))<1e-8 && abs(pointsP1(iN,3)-zTahoe(k))<1e-8)
                    flag = 1;
                    break;
                else
                    flag = 0;
                end
            end
            if(flag)
                break;
            end
        end
        if(flag)
            break;
        end
    end
    if(flag)
        indexFEM = nXY*(k-1) + nxT*(j-1)+ i;
        fprintf(fid,'%d %d\n', iN, indexFEM);
        unstructured2FEM(iN) = indexFEM;
    else
        fprintf(fid,'%d %d\n',iN,-1);
        unstructured2FEM(iN) = -1;
    end
end
fclose(fid);


for k = 1: nzT
    for j = 1: nyT
        for i = 1: nxT
            for iN = 1:nPointsForUnstructured
                if(abs(pointsP1(iN,1)-xTahoe(i))<1e-8 && abs(pointsP1(iN,2)-yTahoe(j))<1e-8 && abs(pointsP1(iN,3)-zTahoe(k))<1e-8)
                    flag = 1;
                    break;
                else
                    flag = 0;
                end
            end
            if(flag)
                indexFEM = nXY*(k-1) + nxT*(j-1)+ i;
                FEM2Unstructured(indexFEM) = iN;
            else
                indexFEM = nXY*(k-1) + nxT*(j-1)+ i;
                FEM2Unstructured(indexFEM) = -1;
            end
        end
    end
end


%% Calculation KBC nodes
nKBC = (nzT)*(nxT);   % bottom all points

% Bottom side
j = 1;
countKBCBot = 0;
for k=1:nzT
    for i=1:nxT
        countKBCBot = countKBCBot + 1;
        indexKBC_Bot(countKBCBot) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

KBC2FEM(1:countKBCBot) = [indexKBC_Bot];

if(nKBC~=length(KBC2FEM))
    disp('nKBC not equal length(KBC2FEM)')
    return;
end

fid = fopen('TahoeKBCnodesBody1.dat','w');
for i = 1 : nKBC
    fprintf(fid,'%d\n', KBC2FEM(i));
end
fclose(fid);

%% Calculation of FBC nodes
nFBC = (nxT)*(nzT); % all points on top

j = nyT;
countFBCT = 0;
for k=1:nzT
    for i=1:nxT
        countFBCT = countFBCT + 1;
        indexFBC_T(countFBCT) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

FBC2FEM(1:countFBCT) = [indexFBC_T];

if(nFBC~=length(FBC2FEM))
    disp('Error: nFBC not equal length(FBC2FEM)')
end

fid = fopen('TahoeFBCnodesBody1.dat','w');
for i = 1 : countFBCT
    fprintf(fid,'%d %d %d %f\n', indexFBC_T(i),1,FEM2Unstructured(indexFBC_T(i)),0.00);
end
fclose(fid);

%% Periodic
% front and back side
nPrBC = (nxT)*(nyT-1) ...  % front side excluding bottom
    + (nxT)*(nyT-1) ... % back side excluding bottom
    + (nzT)*(nyT-1) ... % left excluding bot
    + (nzT)*(nyT-1); % right excluding bot
k = 1;
countPrBCF = 0;
for j=2:nyT
    for i=1:nxT
        countPrBCF = countPrBCF + 1;
        indexPrBC_F(countPrBCF) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

k = nzT;
countPrBCB = 0;
for j=2:nyT
    for i=1:nxT
        countPrBCB = countPrBCB + 1;
        indexPrBC_B(countPrBCB) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

% left right
countKBCL = 0;
i = 1;
for k=1:nzT
    for j=2:nyT
        countKBCL = countKBCL + 1;
        indexKBC_L(countKBCL) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

i = nxT;
countKBCR = 0;
for k=1:nzT
    for j=2:nyT
        countKBCR = countKBCR + 1;
        indexKBC_R(countKBCR) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

PrBC2FEM(1:countPrBCF+countPrBCB+countKBCL+countKBCR) = [indexPrBC_F indexPrBC_B indexKBC_L indexKBC_R ];
if(nPrBC~=length(PrBC2FEM))
    disp('Error: nPrBC not equal length(PrBC2FEM)')
end

if(countPrBCF~=countPrBCB)
    disp('Error: countPrBCF not equal countPrBCB')
end

if(countKBCL~=countKBCR)
    disp('countKBCL not equal countKBCR')
end

fid = fopen('TahoePrBC_FB_nodesBody1.dat','w');
for i = 1 : countPrBCF
    fprintf(fid,'%d %d %d\n',indexPrBC_F(i),FEM2Unstructured(indexPrBC_F(i)),0);
end
for i = 1 : countPrBCB
    fprintf(fid,'%d %d %d\n',indexPrBC_B(i),FEM2Unstructured(indexPrBC_B(i)),1);
end
fclose(fid);

fid = fopen('TahoePrBC_LR_nodesBody1.dat','w');
for i = 1 : countKBCL
    fprintf(fid,'%d %d %d\n',indexKBC_L(i),FEM2Unstructured(indexKBC_L(i)),0);
end
for i = 1 : countKBCR
    fprintf(fid,'%d %d %d\n',indexKBC_R(i),FEM2Unstructured(indexKBC_R(i)),1);
end
fclose(fid);

disp("KBC nodes: " + nKBC)
disp("FBC nodes: " + nFBC)
disp("PrBC-FB nodes: " + (countPrBCF+countPrBCB))
disp("PrBC-LR nodes: " + (countKBCL+countKBCR))
