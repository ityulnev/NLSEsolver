%Sets boundary conditions left and right of radial coordinate E(r)
function [Erf]=set_boundaries(mesh,pulse,medium,Erf,M_fd,mm)


%Calc Phase of wave before boundary e^(i*kr*dr)
switch M_fd.boundary
    case 'open'% rough estimate
        %Left boundary 
        kr_L=log(Erf(3,:)./Erf(2,:))./(1i*mesh.dr).*(mesh.r(1)/mesh.r(2));
        Erf(1,:)=Erf(2,:).*exp(-1i.*kr_L.*mesh.dr);
        kr=log(Erf(end-1,:)./Erf(end-2,:))./(1i*mesh.dr);
        kr=(kr).*mesh.r(end)/mesh.r(end-1); 
%         kr(real(kr)>0)=0+1i.*imag(kr(real(kr)>0));
        Erf(end,:)=Erf(end-1,:).*exp(1i.*(kr).*mesh.dr);%(Eprev(end,:)./Eprev(end-1,:));
    case 'openCorrected'
        %Left
        kr_L=log(Erf(3:4,:)./Erf(2:3,:))./(1i*mesh.dr);
        kr_LNext=2.*kr_L(1,:)-1.*kr_L(2,:);                                 %delKr(2)=delKr(3)-dKr, with dKr=delKr(4)-delKr(3)
        Erf(1,:)=Erf(2,:).*exp(-1i.*kr_LNext.*mesh.dr);
        %Right
        kr=log(Erf(end-2:end-1,:)./Erf(end-3:end-2,:))./(1i*mesh.dr);
        krNext=2.*kr(2,:)-kr(1,:);                                          %delKr(Nr+1)=delKr(Nr)+dKr, with dKr=delKr(Nr)-delKr(Nr-1), dKr is a constant value! b.c. linear function with dr->const
        Erf(end,:)=Erf(end-1,:).*exp(1i.*(krNext).*mesh.dr);                %(Eprev(end,:)./Eprev(end-1,:));
    case 'reflect'
        Erf(1,:)=Erf(2,:);
        Erf(end,:)=Erf(end-1,:);
%         Erf(end,:)=zeros(1,mesh.flength);
    case 'fiber'
        Erf(end,:)=Erf(end-1,:).*((medium.n0^2.*mesh.r(end-1))./(1.45^2.*mesh.r(end)));%Refr index step HCF boundary        
end

end