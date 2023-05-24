% this is for Tahoe - Double layer of triangles

clear all
geometry = 'rectangle';

if(strcmp(geometry,'rectangle'))
    
    if(0)
    delta = 1;
    xs = -36; ys = -1; zs = -4; % left bottom corner
    lengthX = 72;
    lengthY = 1; % THE RESIDUAL IS NAN
    %lengthY = 1.001; % outBItriangle.err
    %lengthY = 1; % outBItriangle.err    
    xe = xs + lengthX;
    ye = ys + lengthY;
    nx = 320;
    ny = 32;
    nPlaneInZ = 36;
    nz = nPlaneInZ;
    leftSide = 0; rightSide = 20; % range of flexbile part - it will not be exactly these points - but it will be very close
    frontSide = -3; backSide = 3;        
    end
    
    if(1)
    delta = 1;
    xs = -0.25; ys = -0.5; zs = -0.25; % left bottom corner
    lengthX = 2*pi+2*0.25;
    lengthY = 0.5; % THE RESIDUAL IS NAN
    %lengthY = 1.001; % outBItriangle.err
    %lengthY = 1; % outBItriangle.err    
    xe = xs + lengthX;
    ye = ys + lengthY;
    nx = 100;
    ny = 20;
    nPlaneInZ = 42;
    nz = nPlaneInZ;
    leftSide = 0; rightSide = 2*pi; % range of flexbile part - it will not be exactly these points - but it will be very close
    frontSide = 0; backSide = pi;    
    end
    
    
    if(0)
    delta = 1;
    xs = -56; ys = -0.25; zs = -5.5; % left bottom corner
    lengthX = 112;
    lengthY = 0.25; % THE RESIDUAL IS NAN
    %lengthY = 1.001; % outBItriangle.err
    %lengthY = 1; % outBItriangle.err    
    xe = xs + lengthX;
    ye = ys + lengthY;
    nx = 500;
    ny = 6;
    nPlaneInZ = 100;
    nz = nPlaneInZ;
    leftSide = -7.5; rightSide = 7.5; % range of flexbile part - it will not be exactly these points - but it will be very close
    frontSide = -3; backSide = 3;    
    end
    
    x = linspace(xs,xe,nx);
    y = linspace(ys,ye,ny);
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    
    % adjusting so that flexible panel exactly starts at leftSide
    dx = (leftSide-xs)/floor((leftSide-xs)/dx);
    nxTahoe = floor((rightSide-leftSide)/dx);
    dx = (rightSide-leftSide)/nxTahoe;
    nxLeft = floor((leftSide-xs)/dx);
    xs = leftSide - nxLeft*dx;
    nxRight = floor((xe-rightSide)/dx);
    xe = rightSide + nxRight*dx;
    x = xs:dx:xe;
    nx = length(x);
    ny = length(y);
    % x-position of the rectangular geom
    newX(1:nx-1) = x(1:end-1);
    newX(end+1:end+ny-1) = x(end)*ones(1,ny-1);
    newX(end+1:end+nx-1) = flip(x(2:end));
    newX(end+1:end+ny-1) = x(1)*ones(1,ny-1);
    
    %y-position of the rectangular geom
    newY(1:nx-1) = y(1)*ones(1,nx-1);
    newY(end+1:end+ny-1) = y(1:end-1);
    newY(end+1:end+nx-1) = y(end)*ones(1,nx-1);
    newY(end+1:end+ny-1) = flip(y(2:end));
    
elseif(strcmp(geometry,'tangentHyperbolicBump'))
    
    k = 0.1;
    delta = k/0.547;
    xs = -55*delta; xe = 88*delta; % left bottom corner
    ys = 0.5;
    lengthY = 0.5;
    lengthX = xe - xs;
    ye = ys + lengthY;
    nx = 2000;
    ny = 30;
    
    leftSide = -1.5; rightSide = 1.5; % range of flexbile part - it will not be exactly these points - but it will be very close
    
    x = linspace(xs,xe,nx);
    y = linspace(ys,ye,ny);
    
    hBump = calculateBumpHeight(k,x);
    
    % x-position of the rectangular geom (counterClockWise)
    newX(1:nx-1) = x(1:end-1);
    newX(end+1:end+ny-1) = x(end)*ones(1,ny-1);
    newX(end+1:end+nx-1) = flip(x(2:end));
    newX(end+1:end+ny-1) = x(1)*ones(1,ny-1);
    
    %y-position of the rectangular geom
    newY(1:nx-1) = y(1)*ones(1,nx-1);
    newY(end+1:end+ny-1) = y(1:end-1);
    newY(end+1:end+nx-1) = y(end)*ones(1,nx-1)+ flip(hBump(2:end));
    newY(end+1:end+ny-1) = flip(y(2:end));
    
end

dx = x(2) - x(1);
dy = y(2) - y(1);
disp(['Boundary layer thickness = ',num2str(delta)]);



%% generating the relation of the nodes for the elements and writing the file unstrucSurfBody1.dat
%       pointsP1(:,1) ----> x coordinate
%       pointsP1(:,2) ----> y coordinate

nP_NewX = length(newX);
zPs = zeros(nP_NewX,1); 
dz = dx;

ze = zs+nPlaneInZ*dz;
% adjusting so that flexible panel exactly starts at leftSide
nzTahoe = floor((backSide-frontSide)/dz);
dz = (backSide-frontSide)/nzTahoe;
nzFront = floor((frontSide-zs)/dz);
zs = frontSide - nzFront*dz;
nzBack = floor((ze-backSide)/dz);
ze = backSide + nzBack*dz;
newTempZ = zs:dz:ze;
ze = newTempZ(end);
nPlaneInZ = length(newTempZ);
disp(['spanwise width will be ',num2str((nPlaneInZ-1)*dz)]);


for i = 1:nPlaneInZ
pointsP1((i-1)*nP_NewX+1:i*nP_NewX,:)=[newX' newY' zPs+newTempZ(i)];
end
newZ = pointsP1(:,3)';
%               plane 1                 plane 2               plane 3

nPoints=length(newX);
nElements=nPoints*2*(nPlaneInZ-1);

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

leftMost = nxLeft+1;
rightMost = nx - nxRight;

z = zs:dz:ze;

frontMost = nzFront + 1;
backMost = nz - nzBack - 1;

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


% Calculation KBC nodes
nKBC = 2*nyT*nzT + 2*(nxT-2)*nyT + (nxT*nzT - 2*nzT - 2*(nxT-2));
i = nxT;
countKBCR = 0;
for k=1:nzT
    for j=1:nyT
        countKBCR = countKBCR + 1;
        indexKBC_R(countKBCR) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

countKBCL = 0;
i = 1;
for k=1:nzT
    for j=1:nyT
        countKBCL = countKBCL + 1;
        indexKBC_L(countKBCL) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

k = 1;
countKBCF = 0;
for j=1:nyT
    for i=2:nxT-1
        countKBCF = countKBCF + 1;
        indexKBC_F(countKBCF) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

k = nzT;
countKBCB = 0;
for j=1:nyT
    for i=2:nxT-1
        countKBCB = countKBCB + 1;
        indexKBC_B(countKBCB) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end
% Bottom side
j = 1;
countKBCBot = 0;
for k=2:nzT-1
    for i=2:nxT-1
        countKBCBot = countKBCBot + 1;
        indexKBC_Bot(countKBCBot) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

KBC2FEM(1:countKBCL+countKBCR+countKBCF+countKBCB+countKBCBot) = [indexKBC_L indexKBC_R indexKBC_F indexKBC_B indexKBC_Bot];

if(nKBC~=length(KBC2FEM))
    disp('nKBC not equal length(KBC2FEM)')
    return;
end

fid = fopen('TahoeKBCnodesBody1.dat','w');
for i = 1 : nKBC
    fprintf(fid,'%d\n', KBC2FEM(i));
end
fclose(fid);


% Calculation of FBC nodes
nFBC = (nxT-2)*(nzT-2);

j = nyT;
countFBCT = 0;
for k=2:nzT-1
    for i=2:nxT-1
        countFBCT = countFBCT + 1;
        indexFBC_T(countFBCT) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end


FBC2FEM(1:countFBCT) = [indexFBC_T];

if(nFBC~=length(FBC2FEM))
        disp('nFBC not equal length(FBC2FEM)')
end

fid = fopen('TahoeFBCnodesBody1.dat','w');
for i = 1 : countFBCT
    fprintf(fid,'%d %d %d %f\n', indexFBC_T(i), 1,FEM2Unstructured(indexFBC_T(i)),0.00);
end
fclose(fid);


disp("KBC nodes: " + nKBC)
disp("FBC nodes: " + nFBC)

%% Functions
%
function [newX,newY] = findSingleLineXY(x,y,smallestDis)
count = 1;
for i = 1:length(x)-1
    distance(count) = sqrt((x(i) - x(i+1))^2 + (y(i) - y(i+1))^2);
    count = count + 1;
end
distance(count) = sqrt((x(length(x)) - x(1))^2 + (y(length(x)) - y(1))^2);
maxDis = max(distance);
if (maxDis>=2*smallestDis) %distance between the points should be large enough
    count=1;
    for i=1:length(x)-1
        % angle for inclined line
        angle= atand((1)*(y(i) - y(i+1))/(x(i) - x(i+1)));
        nNewPoint(i)=floor(distance(i)/smallestDis);% how many grid for line i
        nNewPoint=double(nNewPoint);
        tempDis=double(nNewPoint(i))*smallestDis;
        differ=distance(i)-tempDis;
        if (differ<(smallestDis/2) && differ>0)
            nNewPoint(i)= nNewPoint(i)-1;
        elseif (differ==0)
            nNewPoint(i)= nNewPoint(i)-1;
        end
        newX(count)=x(i);
        newY(count)=y(i);
        countOld=count;
        if (x(i)>x(i+1))%current point is > the next point
            for j=count:count+nNewPoint(i)-1
                count=count+1;
                %so we will minus each distance
                newX(j+1)=newX(j) - smallestDis*cosd(angle);
                % we have to move this smallestDis*cosd(angle) amount
                %in x-direction
            end
        else
            for j=count:count+nNewPoint(i)-1
                newX(j+1)=newX(j) + smallestDis*cosd(angle);
                count=count+1;
            end
        end
        
        if (x(i)<x(i+1))
            for j=countOld:countOld+nNewPoint(i)-1
                newY(j+1)=newY(j) + smallestDis*sind(angle);
            end
        else
            for j=countOld:countOld+nNewPoint(i)-1
                newY(j+1)=newY(j) - smallestDis*sind(angle);
            end
        end
        count=count+1;
    end
end
end

%%
function writeFile(outputFile,totalPoints,nElements)
%fullFileName = fullfile('C:\Users\alsha\Documents\result4', 'unstrucSurfBody1.dat');
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


%%

function h = calculateBumpHeight(k,x)
Sr = 20;
Lr = 0.2;
h = 0.5*k*(tanh(Sr*(x+Lr)) - tanh(Sr*(x-Lr)));
end