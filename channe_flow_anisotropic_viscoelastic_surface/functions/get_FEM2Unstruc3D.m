function FEM2Unstructured = get_FEM2Unstruc3D(pointsP1,xT,yT,zT,eps)

nxT = length(xT);
nyT = length(yT);
nzT = length(zT);
nN = length(pointsP1(:,1));

FEM2Unstructured = zeros(nxT*nyT*nzT,1);

for k = 1: nzT
    for j = 1: nyT
        for i = 1: nxT
            for iN = 1:nN
                if(abs(pointsP1(iN,1)-xT(i))<eps && abs(pointsP1(iN,2)-yT(j))<eps && abs(pointsP1(iN,3)-zT(k))<eps)
                    flag = 1;
                    break;
                else
                    flag = 0;
                end
            end
            if(flag)
                indexFEM = nxT*nyT*(k-1) + nxT*(j-1)+ i;
                FEM2Unstructured(indexFEM) = iN;
            else
                indexFEM = nxT*nyT*(k-1) + nxT*(j-1)+ i;
                FEM2Unstructured(indexFEM) = -1;
            end
        end
    end
end
end