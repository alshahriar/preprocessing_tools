function [unstructured2FEM,HypoMarkerIndex] = get_unstructured2FEM3D_KBCasFBC(unstructured2FEM,pointsP1,xyzHypoKBCRight,indexFEM,eps)
nN = length(pointsP1(:,1));
nNHypo = length(xyzHypoKBCRight(:,1));
HypoMarkerIndex = zeros(nNHypo,1);
count = 0;
maxXUnstruc = max(pointsP1(:,1));
for iN = 1:nN
    for m = 1:nNHypo
        xT = xyzHypoKBCRight(m,1);
        yT = xyzHypoKBCRight(m,2);
        zT = xyzHypoKBCRight(m,3);
        if(maxXUnstruc==pointsP1(iN,1) && abs(pointsP1(iN,2)-yT)<eps && abs(pointsP1(iN,3)-zT)<eps)
            unstructured2FEM(iN) = indexFEM(m);
            count = count + 1;
            HypoMarkerIndex(count) = iN;
            break;
        end
    end
end
disp("Hypo point found: " ...
    + count + " nNHypo: " + nNHypo)
end