function [m1,m2,S,ti,tj,tk] = computeGrads(HGP)

NF = numel(HGP.F)/3;

v1 = HGP.V(HGP.F(:,2),:) - HGP.V(HGP.F(:,1),:) ;
v2 = HGP.V(HGP.F(:,3),:) - HGP.V(HGP.F(:,1),:) ;

Xi = zeros(NF,1);
Yi = Xi;
Yj = Xi;
Xj = sqrt(sum(abs( v1 ).^2,2));
n13 = sqrt(sum(abs( v2 ).^2,2));
n23 = sqrt(sum(abs( HGP.V(HGP.F(:,2),:)-HGP.V(HGP.F(:,3),:) ).^2,2));
Xk = (Xj.^2 + n13.^2 - n23.^2)./(2*Xj);
Yk = sqrt(n13.^2 - Xk.^2);

S = [v1(:,2).*v2(:,3)-v2(:,2).*v1(:,3) -1*(v1(:,1).*v2(:,3)-v2(:,1).*v1(:,3)) v1(:,1).*v2(:,2)-v2(:,1).*v1(:,2)];
S = sqrt(sum(abs( S ).^2,2)) / 2;

temp = 1./(2*S);
m1 = [temp temp temp].*[Yj-Yk Yk-Yi Yi-Yj];
m2 = [temp temp temp].*[Xk-Xj Xi-Xk Xj-Xi];

e_k = Xj+1i*Yj-(Xi+1i*Yi);
e_i = Xk+1i*Yk-(Xj+1i*Yj);
e_j = Xi+1i*Yi-(Xk+1i*Yk);

ti = e_i*1i;
tj = e_j*1i;
tk = e_k*1i;

end

