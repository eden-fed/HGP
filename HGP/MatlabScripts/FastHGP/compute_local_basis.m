function [ e1 e2 ] = compute_local_basis(v1,v2)
%[ e1 e2 ] = compute_local_basis(v1,v2)
%   v1  edge vector 1 per triangle  #F x 3
%   v2  edge vector 2 per triangle  #F x 3
%   e1,e2: two #F x 3  arrays of unit-length vectors, orthogonal to each other; local
  len1 = (sqrt(sum((v1').^2)))';
  e1 = v1./[len1,len1,len1]; 
  v2proj = dot(e1',v2')';
  e2 = v2 - [v2proj,v2proj,v2proj].*e1;
  len2 = (sqrt(sum((e2').^2)))';
  e2 = e2 ./[len2,len2,len2];
end

