%This script set the local gradints and the frames for the *first*
%iteration. In the next iteration we will invoke the setFrames2 script.

HGP.Result.timeVector=[];
HGP.hN = max(max(HGP.halfEdges));
[HGP.Grad.m1, HGP.Grad.m2, HGP.Grad.S, HGP.Grad.ti, HGP.Grad.tj, HGP.Grad.tk] = computeGrads(HGP);

if HGP.calcFramesFromVecField
    HGP.frames = getFramesFromVectorField(HGP);
end