function [tc,dblAc] = compute_perps(V, F, xlb)
% [tc,dblAc] = compute_perps(V, F)
% xlb = index of edge vector to be the x axis direction in the local basis
%  tc = # triangles x 3 array of outward perps to sides (length = side
%  length) in complex form in local coords 
% dblA =  array of doubled triangle areas
  % renaming indices of vertices of triangles for convenience
  i1 = F(1,:); i2 = F(2,:); i3 = F(3,:); 
  % #F x 3 matrices of triangle edge vectors, named after opposite vertices
  v1 = V(i3,:) - V(i2,:);  v2 = V(i1,:) - V(i3,:); v3 = V(i2,:) - V(i1,:);

  vlb1=v1; vlb2=v2;
  if nargin>2
      if xlb==2
          vlb1=v2; vlb2=v3;
      elseif xlb==3
          vlb1=v3; vlb2=v1;
      end
  end
 
  [ e1 e2 ] = compute_local_basis(vlb1,vlb2);
  
  % complex edge vectors in local coordinates
  cv1 = dot(e1',v1')' + 1i*dot(e2',v1')';
  cv2 = dot(e1',v2')' + 1i*dot(e2',v2')';
  cv3 = -cv1-cv2;
  % double triangle area
  % crossprod(a,b) = imag( conj(a)*b), dotprod(a,b)= real(conj(a)*b)
  dblAc = abs(imag( conj(cv1) .* cv2));
  % complex form of the first basis vector
  e1c = e1(:,1) + e1(:,2)*1i;
  
  % clockwise rotated side perps
  % not really needed -- can use cv's directly to compute 
  % the matrix, but leaving for now
  t1r = -1i*cv1;
  t2r = -1i*cv2;
  t3r = -1i*cv3;
  tc = [t1r t2r t3r];
end
