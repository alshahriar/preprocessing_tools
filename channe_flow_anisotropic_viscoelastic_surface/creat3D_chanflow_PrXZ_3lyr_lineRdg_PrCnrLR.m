% this is for Tahoe - Double layer of triangles

clear all
geometry = 'rectangle';

if(strcmp(geometry,'rectangle'))
    
    if(1)
        delta = 1;
        xOffset = 0.25;
        zOffset = 0.25;
        xs = -xOffset; ys = -0.5; zs = -zOffset; % left bottom corner
        lengthX = 2*pi;

        %layers
        l1 = 0.2;
        l2 = 0.1;
        l3 = 0.2;
        nn1 = 8;
        nn2 = 4;
        nn3 = 8;
        nRidgeZ = 5;
        nRidgeX = 20;
        lengthY = l1+l2+l3;    
        lengthZ = pi;
        %lengthY = 1.001; % outBItriangle.err
        %lengthY = 1; % outBItriangle.err    
        xe = xs + lengthX + 2*xOffset;
        ye = ys + lengthY;
        nx = 130;
        leftSide = 0; rightSide = 2*pi; % range of flexbile part - it will not be exactly these points - but it will be very close
        frontSide = 0; backSide = pi;    
    end
    
    x = linspace(xs,xe,nx);
    y1  = linspace(ys,ys+l1,nn1);
    y2 = linspace(ys+l1,ys+l1+l2,nn2);
    y3 =  linspace(ys+l1+l2,ys+l1+l2+l3,nn3);
    y = [y1 y2(2:nn2) y3(2:nn3)];
    dx = x(2) - x(1);
    
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
disp(['Boundary layer thickness = ',num2str(delta)]);



%% generating the relation of the nodes for the elements and writing the file unstrucSurfBody1.dat
%       pointsP1(:,1) ----> x coordinate
%       pointsP1(:,2) ----> y coordinate

nP_NewX = length(newX);
zPs = zeros(nP_NewX,1); 
dz = dx;

ze = zs + lengthZ + 2*zOffset;
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

z = zs:dz:ze;
nz = length(z);

for k=1:nz
    if(z(k)>frontSide && z(k)<(frontSide+dz))
        frontMost=k; %index of P2 point
    elseif (abs(z(k)-frontSide)<1e-8)
        frontMost=k;break;
    end
end

for k=1:nz
    if(z(k)<backSide && z(k)>(backSide-dz))
        backMost=k;% index of another point
    elseif(abs(z(k)-backSide)<1e-8)
        backMost=k;break;
    end
end

% frontMost = nzFront + 1;
% backMost = nz - nzBack - 1;

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

if(nFBC*3>20000)
    disp("nFBC*3< 20000 - hard coding restriction in TAHOE SRC")
end

j = nyT;
countFBCT = 0;
for k=1:nzT
    for i=1:nxT
        countFBCT = countFBCT + 1;
        indexFBC_T(countFBCT) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

FBC2FEM(:) = [indexFBC_T];

if(nFBC~=length(FBC2FEM))
        disp('Error: nFBC not equal length(FBC2FEM)')
end

fid = fopen('TahoeFBCnodesBody1.dat','w');
for i = 1 : countFBCT
    fprintf(fid,'%d %d %d %f\n', indexFBC_T(i),1,FEM2Unstructured(indexFBC_T(i)),0.00);
end
fclose(fid);

% figure;
% hold on
fid = fopen('FBCxyzLocFromMatlab.dat','w');
fprintf(fid,'%d\n',countFBCT); 
for i = 1 : countFBCT
    fprintf(fid,'%15e %15e %15e\n',pointTemp(indexFBC_T(i),1),pointTemp(indexFBC_T(i),2),pointTemp(indexFBC_T(i),3)); 
    %plot3(pointTemp(indexFBC_T(i),1),pointTemp(indexFBC_T(i),2),pointTemp(indexFBC_T(i),3),'.')
end
fclose(fid);

%% Periodic 
% front and back side
nPrBC = (nxT)*(nyT-1) ...  % front side excluding bottom
      + (nxT)*(nyT-1) ... % back side excluding bottom
      + (nzT-2)*(nyT-1) ... % left excluding bot left right
      + (nzT-2)*(nyT-1); % right excluding bot left right
  
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
for k=2:nzT-1
    for j=2:nyT
        countKBCL = countKBCL + 1;
        indexKBC_L(countKBCL) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

i = nxT;
countKBCR = 0;
for k=2:nzT-1
    for j=2:nyT
        countKBCR = countKBCR + 1;
        indexKBC_R(countKBCR) = i + (j-1)*nxT + (k-1)*nxT*nyT;
    end
end

PrBC2FEM(:) = [indexPrBC_F indexPrBC_B indexKBC_L indexKBC_R ];
if(nPrBC~=length(PrBC2FEM))
        disp('nPrBC not equal length(PrBC2FEM)')
end

if(countPrBCF~=countPrBCB)
        disp('Error: countPrBCF not equal countPrBCB')
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

disp("FBC KBC nodes: " + nFBC + " " + nKBC)

disp("PrBC-LR nodes: PrBCnodesX: " + (countKBCL+countKBCR))
disp("PrBC-FB nodes: PrBCnodesZ: " + (countPrBCF+countPrBCB))
%% Element-wise layers


% nXY = nxT*nyT;
% Ee = (nzT-1)*(nyT-1)*(nxT-1);
% elmt = zeros(Ee,8);
% ie = 0;
% for k = 1:(nzT-1)
%     for j = 1:(nyT-1)
%         for i = 1:(nxT-1)
%             ie = ie + 1;
%             elmt(ie,:) = [nXY*(k-1) + nxT*(j-1)+ i      nXY*(k-1) + nxT*(j-1)+i+1    nXY*(k-1) + nxT*j + i+1   nXY*(k-1) + nxT*j+ i ...
%                           nXY*(k-0) + nxT*(j-1)+ i      nXY*(k-0) + nxT*(j-1)+i+1    nXY*(k-0) + nxT*j + i+1   nXY*(k-0) + nxT*j+ i ];            
% %            elmt(ie,:) = [nxT*(j-1)+i     nxT*(j-1)+i+1    nxT*j+i+1   nxT*j+i];
%         end
%     end
% end

init_nzDr = 3;
nzDR = init_nzDr;
init_gapZ_DR = 9;
gapZ_DR = init_gapZ_DR;

minNumberOfRidges = 4;
changeGapFlag = 1;
changeSizeFlag = 0;


nxElm = nxT-1;
nyElm = nyT-1;
nzElm = nzT-1;
nnELyr1 = (nn1-1);
nnELyr2 = (nn2-1);
nnELyr3 = (nn3-1);
nXYElm = nxElm*nyElm;
elmtCheck = -1*ones(Ee,1);

% Dermal Ridges on the 1st layer
nxDR = nxElm; % size of each Dermal Ridge
nyDR = nnELyr1;

NnXDR = 1; % How many Dermal Ridges in x dir
gapX_DR = 0;
x_DR_S = 0;
z_DR_S = 0;
            
min_NnZDR = 4;
min_nzDR = 4;
min_gapZ_DR = min_nzDR+1;
max_NnZDR = floor((nzElm-min_nzDR)/min_gapZ_DR);
max_gapZ_DR = floor(nzElm/min_gapZ_DR);
max_nzDR = floor(nzElm/min_NnZDR - 1);

maxEmltZ = z_DR_S + (min_NnZDR-1)*min_gapZ_DR + max_nzDR;
if(maxEmltZ>nzElm)
    disp("failed:  max_nzDR")
end

maxEmltZ = z_DR_S + (max_NnZDR-1)*min_gapZ_DR + min_nzDR;
if(maxEmltZ>nzElm)
    disp("failed:  max_NnZDR")
end

maxEmltZ = z_DR_S + (min_NnZDR-1)*max_gapZ_DR + min_nzDR;
if(maxEmltZ>nzElm)
    disp("failed:  max_gapZ_DR")
end

% finding all the combinations that give symmetric configuration
clear set
count = 0;
for NnZDR = min_NnZDR: max_NnZDR
    for gapZ_DR = min_gapZ_DR:max_gapZ_DR
        for nzDR = min_nzDR:gapZ_DR-1
            nzElmOccupied = NnZDR*gapZ_DR;
            if((nzElm - nzElmOccupied) == 0)
                count = count + 1;
                disp(count + " : Nn of DR slices = " + NnZDR + "; gapZ_DR = " + gapZ_DR + "; nzDR = " + nzDR + ";")
                set(count,1:3) = [NnZDR gapZ_DR nzDR];
            end            
        end
    end
end

% selecting final values of the configurations 
userSelect = input(['which one? 1-',num2str(count),' :']);
[NnZDR] = set(userSelect,1);
[gapZ_DR] = set(userSelect,2);
[nzDR] = set(userSelect,3);


% checking if occupied elements are more than available elemetns
maxEmltX = x_DR_S + (NnXDR-1)*gapX_DR + nxDR;
if(maxEmltX>nxElm)
    disp("too many DR in x-dir")
    quit;
end

maxEmltZ = z_DR_S + (NnZDR-1)*gapZ_DR + nzDR;
if(maxEmltZ>nzElm)
    disp("too many DR in z-dir")
    quit;
end

nElmtEachDR = nxDR*nyDR*nzDR;
totalElmtForDR1 = NnXDR*NnZDR*nElmtEachDR;
elmtListDR = zeros(nElmtEachDR,1);
DR_ElmtXY1 = zeros(NnXDR,NnZDR,nElmtEachDR);
DR_ElmtList1 = zeros(totalElmtForDR1,1);


total_count  = 0;
for i = 1: NnXDR
    for k = 1:NnZDR
        ie = 0;
        for kk = 1:nzDR
            for jj = 1:nyDR
                for ii = 1:nxDR         
                    iiF = x_DR_S + (i-1)*gapX_DR + ii;
                    jjF = jj;
                    kkF = z_DR_S + (k-1)*gapZ_DR + kk;
                    elmtID = nXYElm*(kkF-1) + nxElm*(jjF-1)+ iiF;                    
                    if(elmtCheck(elmtID) == -1)
                        ie = ie + 1;
                        elmtListDR(ie) = elmtID;
                        elmtCheck(elmtID) = 4;
                    end
                end
            end
        end
        total_count = total_count + ie;
        DR_ElmtXY1(i,k,:) = elmtListDR;     
    end
end

if(total_count ~= totalElmtForDR1)
    disp("Error: DR 1 have problems")
end

iCount = 0;
for i = 1: NnXDR
    for k = 1:NnZDR
        for e = 1:nElmtEachDR
            iCount = iCount + 1;
            DR_ElmtList1(iCount) = DR_ElmtXY1(i,k,e);
        end
    end
end


% DR on the 2nd Layer
nyDR = nnELyr2;
nElmtEachDR = nxDR*nyDR*nzDR;
totalElmtForDR2 = NnXDR*NnZDR*nElmtEachDR;
elmtListDR = zeros(nElmtEachDR,1);
DR_ElmtXY2 = zeros(NnXDR,NnZDR,nElmtEachDR);
DR_ElmtList2 = zeros(totalElmtForDR2,1);

total_count = 0;
for i = 1: NnXDR
    for k = 1:NnZDR
        ie = 0;
        for kk = 1:nzDR
            for jj = 1:nyDR
                for ii = 1:nxDR
                    iiF = x_DR_S + (i-1)*gapX_DR + ii;
                    jjF = nnELyr1+jj;
                    kkF = z_DR_S + (k-1)*gapZ_DR + kk;
                    elmtID = nXYElm*(kkF-1) + nxElm*(jjF-1)+ iiF;
                    if(elmtCheck(elmtID) == -1)
                        ie = ie + 1;
                        elmtListDR(ie) = elmtID;
                        elmtCheck(elmtID) = 5;
                    end
                end
            end
        end
        total_count = total_count + ie;
        DR_ElmtXY2(i,k,:) = elmtListDR;
    end
end

if(total_count ~= totalElmtForDR2)
    disp("Error: DR 2 have problems")
end

iCount = 0;
for i = 1: NnXDR
    for k = 1:NnZDR
        for e = 1:nElmtEachDR
            iCount = iCount + 1;
            DR_ElmtList2(iCount) = DR_ElmtXY2(i,k,e);
        end
    end
end

% 3 Layers
% Layer 1
nEL1 = nxElm*nzElm*nnELyr1 - totalElmtForDR1;
elmtGrp1 = zeros(nEL1,1);
ie = 0;
for k = 1:nzElm
    for j = 1:nnELyr1
        for i = 1:nxElm
            elmtID = nXYElm*(k-1) + nxElm*(j-1)+ i;
            if(elmtCheck(elmtID) == -1)
                ie = ie + 1;
                elmtGrp1(ie) = elmtID;
                elmtCheck(elmtID) = 1;
            end
        end
    end
end
if(ie~=nEL1)
    disp("Error: grp 1")
end
% Layer 2
nEL2 = nxElm*nzElm*nnELyr2 - totalElmtForDR2;
elmtGrp2 = zeros(nEL2,1);
ie = 0;
for k = 1:nzElm
    for j = nnELyr1+1:nnELyr1+nnELyr2
        for i = 1:nxElm
            elmtID =  nXYElm*(k-1) + nxElm*(j-1)+ i;            
            if(elmtCheck(elmtID) == -1)
                ie = ie + 1;
                elmtGrp2(ie) = elmtID;
                elmtCheck(elmtID) = 2;
            end
        end
    end
end
if(ie~=nEL2)
    disp("Error: grp 2")
end

% Layer 3
nEL3 = nxElm*nzElm*nnELyr3;
elmtGrp3 = zeros(nEL3,1);
ie = 0;
for k = 1:nzElm
    for j = nnELyr1+nnELyr2+1:nnELyr1+nnELyr2+nnELyr3
        for i = 1:nxElm
            elmtID = nXYElm*(k-1) + nxElm*(j-1)+ i;
            if(elmtCheck(elmtID) == -1)
                ie = ie + 1;
                elmtGrp3(ie) = elmtID;
                elmtCheck(elmtID) = 3;
            end
        end
    end
end
if(ie~=nEL3)
    disp("Error: grp 3")
end

if(sum(elmtCheck==-1)>0)
    disp("Error: some elements has not been assigned")
end


%% Saving the list of elements

fid = fopen('Tahoe_Group1_ElemList_Body1.dat','w');
fprintf(fid,'%d\n',nEL1);
for i = 1 : nEL1
    fprintf(fid,'%d %d\n',elmtGrp1(i),elmtCheck(elmtGrp1(i)));
end
fclose(fid);

fid = fopen('Tahoe_Group2_ElemList_Body1.dat','w');
fprintf(fid,'%d\n',nEL2);
for i = 1 : nEL2
    fprintf(fid,'%d %d\n',elmtGrp2(i),elmtCheck(elmtGrp2(i)));
end
fclose(fid);

fid = fopen('Tahoe_Group3_ElemList_Body1.dat','w');
fprintf(fid,'%d\n',nEL3);
for i = 1 : nEL3
    fprintf(fid,'%d %d\n',elmtGrp3(i),elmtCheck(elmtGrp3(i)));
end
fclose(fid);

fid = fopen('Tahoe_Group4_ElemList_Body1.dat','w');
fprintf(fid,'%d\n',totalElmtForDR1);
for i = 1 : totalElmtForDR1
    fprintf(fid,'%d %d\n',DR_ElmtList1(i),elmtCheck(DR_ElmtList1(i)));
end
fclose(fid);

fid = fopen('Tahoe_Group5_ElemList_Body1.dat','w');
fprintf(fid,'%d\n',totalElmtForDR2);
for i = 1 : totalElmtForDR2
    fprintf(fid,'%d %d\n',DR_ElmtList2(i),elmtCheck(DR_ElmtList2(i)));
end
fclose(fid);

%% CCNS Way of writing

if(Ee ~= nEL1+nEL2+nEL3+totalElmtForDR1+totalElmtForDR2)
    disp("something is wrong with total number of Elmt")
end

disp("Info needed for tahoe");
disp(nEL1 + " " + nEL2 + " " + nEL3 + " " + totalElmtForDR1 + " " + totalElmtForDR2)
disp("Nn Ee: " + Nn + "  " + Ee)
disp("nxTahoe nzTahoe: " + nxT + "  " + nzT)
disp("unStructured Nn Ee: " + nPlaneInZ*nPoints + "  " + nElements)
disp("strideX: " + lengthX + " - strideZ: " + lengthZ)
disp("xTahoe: " + xTahoe(1) + "-" + xTahoe(end))
disp("zTahoe: " + zTahoe(1) + "-" + zTahoe(end))

fid = fopen('TahoeEltGroupsBody1.dat','w');
for i = 1 : nEL1
    ie = elmtGrp1(i);
    fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',i,elmt(ie,1),elmt(ie,2),elmt(ie,3),elmt(ie,4),elmt(ie,5),elmt(ie,6),elmt(ie,7),elmt(ie,8));
end
for i = 1 : nEL2
    ie = elmtGrp2(i);
    fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',i,elmt(ie,1),elmt(ie,2),elmt(ie,3),elmt(ie,4),elmt(ie,5),elmt(ie,6),elmt(ie,7),elmt(ie,8));
end
for i = 1 : nEL3
    ie = elmtGrp3(i);
    fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',i,elmt(ie,1),elmt(ie,2),elmt(ie,3),elmt(ie,4),elmt(ie,5),elmt(ie,6),elmt(ie,7),elmt(ie,8));
end
for i = 1 : totalElmtForDR1
    ie = DR_ElmtList1(i);
    fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',i,elmt(ie,1),elmt(ie,2),elmt(ie,3),elmt(ie,4),elmt(ie,5),elmt(ie,6),elmt(ie,7),elmt(ie,8));
end
for i = 1 : totalElmtForDR2
    ie = DR_ElmtList2(i);
    fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',i,elmt(ie,1),elmt(ie,2),elmt(ie,3),elmt(ie,4),elmt(ie,5),elmt(ie,6),elmt(ie,7),elmt(ie,8));
end
fclose(fid);




%%
filename = 'skinSkelaton.gif';

Plot_Mesh3D(pointTemp,elmt(elmtGrp3,:),'g');
ylim([-0.5 0])
view(-40,20)
hold on
frame = getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',inf);

Plot_Mesh3D(pointTemp,elmt(DR_ElmtList2,:),'m');
frame = getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif','WriteMode','append');

Plot_Mesh3D(pointTemp,elmt(elmtGrp2,:),'b');
frame = getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif','WriteMode','append');

Plot_Mesh3D(pointTemp,elmt(DR_ElmtList1,:),'y');
frame = getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif','WriteMode','append');

Plot_Mesh3D(pointTemp,elmt(elmtGrp1,:),'r');
frame = getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif','WriteMode','append');

%% Functions
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
