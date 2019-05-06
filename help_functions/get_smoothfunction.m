%Take matrix (r,f) and calculate smoothing functions for Maximum in r and f
%A smoothfunction is used to get rid of sharp edges that cuase problems in
%fourier transform e.g. FFT(smoothfct*stepfct)
function [r_smoothfct,f_smoothfct]=get_smoothfunction(mesh,Erf,type)
%find max
[maxf,imaxf]=max(max(abs(Erf).^2,[],1),[],2);
[maxr,imaxr]=max(max(abs(Erf).^2,[],2),[],1);
%index for values below threshold in r and f
rbound=find(abs(Erf(:,imaxf).^2)<maxr.*1e-9,1);
fbound=find(abs(Erf(imaxr,:).^2)>maxf.*1e-9,1);

halfwidthR=mesh.dr*abs(rbound)+(mesh.rmin-mesh.dr);
halfwidthF=mesh.df*abs(mesh.fbound-fbound);

if isempty(rbound)
    halfwidthR=inf;
end
if isempty(fbound)
    halfwidthF=0;
end
%calculate the smoothing function
switch type
    case 'supergaussian'
    r_smoothfct=calc_supergaussian(mesh.r,halfwidthR,10);
    r_smoothfct=r_smoothfct';
    f_smoothfct=calc_supergaussian(mesh.f,halfwidthF,10);
end

end