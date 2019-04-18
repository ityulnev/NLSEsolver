%Sets boundary conditions left and right of radial coordinate E(r)
function [Erf]=set_boundaries(mesh,medium,Erf,M_fd,Eprev,mm)

%Left boundary
Erf(1,:)=Erf(2,:);
kr_L=log(Eprev(1,:)./Eprev(2,:))./(1i*mesh.dr);
% kr_L=kr_L.*(mesh.r(1)/mesh.r(2));
% kr_L(real(kr_L)<0)=0+1i.*imag(kr_L(real(kr_L)<0));
Erf(1,:)=Erf(2,:).*exp(1i.*kr_L.*mesh.dr);
%Calc Phase of reflected wave e(i*kr*dr)
kr=log(Erf(end-1,:)./Erf(end-2,:))./(1i*mesh.dr);


kr=(kr).*mesh.r(end)/mesh.r(end-1); 
% kr(real(kr)<0)=0+1i.*imag(kr(real(kr)<0));
%Right boundary
switch M_fd.boundary
    case 'openR'
        Erf(end,:)=Erf(end-1,:).*(Eprev(end,:)./Eprev(end-1,:)).*(mesh.r(end-1)/mesh.r(end));
    case 'open'
        Erf(end,:)=Erf(end-1,:).*exp(1i.*(kr).*mesh.dr);%(Eprev(end,:)./Eprev(end-1,:));

    case 'reflect'
        Erf(end,:)=zeros(1,mesh.flength);
    case 'fiber'
        Erf(end,:)=Erf(end-1,:).*((medium.n0^2.*mesh.r(end-1))./(1.45^2.*mesh.r(end)));%Refr index step HCF boundary
        
end

end