
% create T_h - harmonic basis matrix for each triangle : dim: 3 X n X |F|
HarmonicBasisR1=full(HarmonicBasis(FastHGP.F(:,1),:)); % first row of T_h in each triangle dim: |F| X n 
HarmonicBasisR2=full(HarmonicBasis(FastHGP.F(:,2),:)); % second row of T_h in each triangle dim: |F| X n
HarmonicBasisR3=full(HarmonicBasis(FastHGP.F(:,3),:)); % third row of T_h in each triangle dim: |F| X n

% create T_fz and T_fbz - get the f_z of each triangle : dim: 1 X 3 X |F| (we use |F| X 3)

[tc,dblAc] = compute_perps(FastHGP.V, FastHGP.F', 3); %tc = |F| x 3 array in complex form, dblA =  array of doubled triangle areas
cf = -(1.0/2.0)./ dblAc;
t1r = tc(:,1); t2r =  tc(:,2); t3r =  tc(:,3);

J_fz=[cf.*conj(t1r) cf.*conj(t2r) cf.*conj(t3r)];% dim:|F| x 3
J_fbz=[cf.*t1r cf.*t2r cf.*t3r];% dim:|F| x 3

% create T_h_fz and T_h_fbz dim: 1 X n X |F| (we use |F| X n)
FastHGP.J_fz=bsxfun(@times, J_fz(:,1), HarmonicBasisR1)+bsxfun(@times, J_fz(:,2), HarmonicBasisR2)+bsxfun(@times, J_fz(:,3), HarmonicBasisR3);
FastHGP.J_fbz=bsxfun(@times, J_fbz(:,1), HarmonicBasisR1)+bsxfun(@times, J_fbz(:,2), HarmonicBasisR2)+bsxfun(@times, J_fbz(:,3), HarmonicBasisR3);

%get array of normalized triangles area
FastHGP.Area=dblAc/2;
FastHGP.Area=FastHGP.Area./sum(FastHGP.Area);

clear HarmonicBasis HarmonicBasisR1 HarmonicBasisR2 HarmonicBasisR3 tc dblAc cf t1r t2r t3r J_fz J_fbz
