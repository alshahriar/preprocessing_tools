function unstructured2FEM = get_unstructured2FEM3D(pointsP1,xT,yT,zT,eps)

nxT = length(xT);
nyT = length(yT);
nzT = length(zT);
nN = length(pointsP1(:,1));

unstructured2FEM = zeros(nN,1);

for iN = 1:nN
    for k = 1:nzT
        for j = 1: nyT
            for i = 1: nxT
                if(abs(pointsP1(iN,1)-xT(i))<eps && abs(pointsP1(iN,2)-yT(j))<eps && abs(pointsP1(iN,3)-zT(k))<eps)
                    flag = 1; % this point is on the edge
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
        indexFEM = nxT*nyT*(k-1) + nxT*(j-1)+ i;
        unstructured2FEM(iN) = indexFEM;
    else
        unstructured2FEM(iN) = -1;
    end
end
end