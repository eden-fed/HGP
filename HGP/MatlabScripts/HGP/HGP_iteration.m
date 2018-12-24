HGP.itTime = tic;
cvx_status = 'null';

if isequal(HGP.W,sparse(length(HGP.W),length(HGP.W)))
    HGP.internalConstraints = [];
else
    HGP.internalConstraints = HGP.W(HGP.rowsToTake,:);
end

HGP.seamConstraintsX = HGP.WSeamX(1:HGP.rNum,:);
HGP.seamConstraintsY = HGP.WSeamY(1:HGP.rNum,:);
HGP.wu1 = HGP.WSeamU1(1:HGP.rNum,:);
HGP.wu2 = HGP.WSeamU2(1:HGP.rNum,:);
HGP.wv1 = HGP.WSeamV1(1:HGP.rNum,:);
HGP.wv2 = HGP.WSeamV2(1:HGP.rNum,:);

[HGP.nrow, HGP.ncol] = size(HGP.internalConstraints);
HGP.harmonicInternal = [HGP.internalConstraints sparse(HGP.nrow,HGP.ncol); sparse(HGP.nrow,HGP.ncol) HGP.internalConstraints];
HGP.seamW = [HGP.seamConstraintsX+HGP.wu1 HGP.wv1; HGP.wu2 HGP.seamConstraintsY+HGP.wv2];

HGP.HarmonicMat2 =[HGP.harmonicInternal;HGP.seamW];

cvx_begin 
cvx_solver_settings( 'MSK_IPAR_PRESOLVE_USE', 'MSK_PRESOLVE_MODE_OFF' )
variable UV(HGP.hN,2);

t1=HGP.Grad.m1(:,1).*UV(HGP.halfEdges(:,1),1)+HGP.Grad.m1(:,2).*UV(HGP.halfEdges(:,2),1)+HGP.Grad.m1(:,3).*UV(HGP.halfEdges(:,3),1);
t2=HGP.Grad.m2(:,1).*UV(HGP.halfEdges(:,1),2)+HGP.Grad.m2(:,2).*UV(HGP.halfEdges(:,2),2)+HGP.Grad.m2(:,3).*UV(HGP.halfEdges(:,3),2);
t3=HGP.Grad.m1(:,1).*UV(HGP.halfEdges(:,1),2)+HGP.Grad.m1(:,2).*UV(HGP.halfEdges(:,2),2)+HGP.Grad.m1(:,3).*UV(HGP.halfEdges(:,3),2);
t4=HGP.Grad.m2(:,1).*UV(HGP.halfEdges(:,1),1)+HGP.Grad.m2(:,2).*UV(HGP.halfEdges(:,2),1)+HGP.Grad.m2(:,3).*UV(HGP.halfEdges(:,3),1);

fz_barReal = 0.5 * (t1-t2);
fz_barImag = 0.5 * (t3+t4);
fzReal = 0.5 * ( t2 + t1 );
fzImag = 0.5 * ( t3 - t4 );

fz_bar = fz_barReal+1i*fz_barImag;
fz = fzReal+1i*fzImag;

E_SOFT = norm(HGP.HarmonicMat2*[UV(:,1);UV(:,2)]);

minimize E_SOFT
subject to

%Rotation constraints
HGP.rotMatrix*[UV(:,1);UV(:,2)] == 0;
%Lipman 12 constraints
real( fz(HGP.BVFaces).*HGP.frames(HGP.BVFaces) ) >= abs(fz_bar(HGP.BVFaces))+0.01; 

%Fix one cone to the origin
UV(HGP.coneIndices(1),1) == 0;
UV(HGP.coneIndices(1),2) == 0;
cvx_end
HGP.fz = fz;
HGP.UV=UV;
clear t1 t2 t3 t4 fz_barReal fz_barImag fzReal fzImag fz_bar fz UV;
clear cvx_*;
HGP.itTime = toc(HGP.itTime);
HGP.Result.timeVector = [HGP.Result.timeVector HGP.itTime];