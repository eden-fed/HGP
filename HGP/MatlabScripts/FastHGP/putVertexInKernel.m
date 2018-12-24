function [ newUV ] = putVertexInKernel( UV, oneRing )
% UV - 2 X 1 array of UV values of the vertex
% oneRing - 2 X size of neighbors of the vertex - ordered counterclockwise

crossSign = @(v,u) (v(1,:).*u(2,:)-u(1,:).*v(2,:));

edges=circshift(oneRing,-1, 2)-oneRing;
oneRingToUV=UV-oneRing;
crossEps=min(abs(crossSign(edges,oneRingToUV)))/100;

cvx_begin quiet
variable newUV(2,1);

oneRingToUV1=newUV(1)-oneRing(1,:);
oneRingToUV2=newUV(2)-oneRing(2,:);
oneRingToUV=[oneRingToUV1;oneRingToUV2];

minimize norm(newUV-UV)
subject to

crossSign(edges,oneRingToUV)>=crossEps;

cvx_end

if(~all(crossSign(edges,oneRingToUV)>=0))%if CVX falied, keep original UV
    newUV=UV;
end

end

