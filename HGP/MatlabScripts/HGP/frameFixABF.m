function [ newOneRingAngles ] = frameFixABF( HGP )
%calculate new angles to the "bad" one ring
FOneRing = HGP.F(HGP.FrameFix.badOneRing,:);
[nr,nc]=size(FOneRing);
cvx_begin
variable newOneRingAngles(nr,nc);
E_ABF = norm( ((newOneRingAngles-HGP.FrameFix.originalOneRingAngle)./HGP.FrameFix.originalOneRingAngle)' , 'fro' ) ;

minimize E_ABF
subject to
newOneRingAngles >= 0.0001;
sum(newOneRingAngles(:,1)) == HGP.FrameFix.oneRingConeAngle;

newOneRingAngles(:,1) + newOneRingAngles(:,2) + newOneRingAngles(:,3) == pi*ones(nr,1);
cvx_end

newOneRingAngles = newOneRingAngles(:,1);
end

