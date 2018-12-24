function [frames] = getFramesFromVectorField(HGP)
% Calculate the frames from a given vector field
% HGP is a struct that contains the faces,vertices and vf of the mesh.

e1 = HGP.V(HGP.F(:,2),:) - HGP.V(HGP.F(:,1),:);
e1 = e1./repmat(sqrt(sum(abs( e1 ).^2,2)),1,3);
   
e2 = HGP.V(HGP.F(:,3),:) - HGP.V(HGP.F(:,1),:);
e2 = e2./repmat(sqrt(sum(abs( e2 ).^2,2)),1,3);
   
n = cross(e1,e2);
n = n./repmat(sqrt(sum(abs( n ).^2,2)),1,3);
   
vec1 = cross(e1,HGP.vectorFieldKv1);
normVec1 = sqrt(sum(abs( vec1 ).^2,2));
normVec1(normVec1 == 0) = 1;
vec1 = vec1./repmat(normVec1,1,3);
signSin = vec1(:,1).*n(:,1) + vec1(:,2).*n(:,2) + vec1(:,3).*n(:,3);
   
cosAng = HGP.vectorFieldKv1(:,1).*e1(:,1) + HGP.vectorFieldKv1(:,2).*e1(:,2) + HGP.vectorFieldKv1(:,3).*e1(:,3); 
sinAng = sqrt(abs(1-cosAng.*cosAng));
frames = cosAng + 1i*(signSin.*sinAng);

end

