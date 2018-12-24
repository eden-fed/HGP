%This script set the frames for the next HGP iteration
%

HGP.frames = abs(HGP.fz)./HGP.fz;
if HGP.FrameFix.status == 0
    HGP.FrameFix.fzz = (1./(4*HGP.Grad.S(HGP.FrameFix.badOneRing))).*(conj(HGP.Grad.ti(HGP.FrameFix.badOneRing)).*HGP.FrameFix.embeddedOneRing(:,1) + conj(HGP.Grad.tj(HGP.FrameFix.badOneRing)).*HGP.FrameFix.embeddedOneRing(:,2) + conj(HGP.Grad.tk(HGP.FrameFix.badOneRing)).*HGP.FrameFix.embeddedOneRing(:,3));
    HGP.FrameFix.rr = abs(HGP.fz(HGP.FrameFix.badOneRing(1)))./HGP.frames(HGP.FrameFix.badOneRing(1))/HGP.FrameFix.fzz(1);
    HGP.FrameFix.localFrames = HGP.FrameFix.fzz.*HGP.FrameFix.rr.*HGP.FrameFix.rotationOffset;
    HGP.FrameFix.localFrames = abs(HGP.FrameFix.localFrames)./HGP.FrameFix.localFrames;
    HGP.frames(HGP.FrameFix.badOneRing)=HGP.FrameFix.localFrames;    
end








