%Sets boundary conditions left and right of radial coordinate E(r)
%default = "open" is linear extrapolation to r(1) and r(end)
function [Erf]=set_boundaries(mesh,pulse,medium,Erf,M_fd,mm,Eprev)
%% Left Boundary
switch M_fd.boundary{1,1}
    case 'constDiv'
        Erf(1,:)=Erf(1,:);
    case 'const'
        Erf(1,:)=Erf(2,:);
    case 'open'
        [~,fleft_smoothfct]=get_smoothfunction(mesh,Erf(2,:),'supergaussian');
        dE=log(Erf(3:4,:)./Erf(2:3,:))./(1i*mesh.dr);
        dE(isnan(dE))=0;
        dE(isinf(dE))=0;
        kr_LNext=2.*dE(1,:)-1.*dE(2,:);
        kr_LNext=kr_LNext.*fleft_smoothfct;
        Erf(1,:)=Erf(2,:).*exp(-1i.*kr_LNext.*mesh.dr);
    case 'open2'
        [~,fleft_smoothfct]=get_smoothfunction(mesh,Erf(2,:),'supergaussian');
        dE=Erf(3:4,:)./Erf(2:3,:);
        dE(isnan(dE))=0;
        dE(isinf(dE))=0;
        dENext=2.*dE(1,:)-1.*dE(2,:);
        dENext=dENext.*fleft_smoothfct;         
        Erf(1,:)=Erf(2,:).*dENext;  
end
%% Right Boundary
switch M_fd.boundary{1,2}
    case 'const'
        Erf(end,:)=Erf(end,:);
    case 'open'
        dkr=log(Erf(end-2:end-1,:)./Erf(end-3:end-2,:))./(1i*mesh.dr);
        dkr(isnan(dkr))=0;
        dkr(isinf(dkr))=0;
        krNext=2.*dkr(2,:)-dkr(1,:);                                          %delKr(Nr+1)=delKr(Nr)+dKr, with dKr=delKr(Nr)-delKr(Nr-1), dKr is a constant value! b.c. linear function with dr->const
        [~,fright_smoothfct]=get_smoothfunction(mesh,Erf(end,:),'supergaussian');
        krNext=krNext.*fright_smoothfct;
        Erf(end,:)=Erf(end-1,:).*exp(1i.*(krNext).*mesh.dr);                      
    case 'open2'
        dkr2=(unwrap(angle(Erf(end-2,:)))-unwrap(angle(Erf(end-3,:))))./mesh.dr;
        dkr2=[dkr2;(unwrap(angle(Erf(end-1,:)))-unwrap(angle(Erf(end-2,:))))./mesh.dr];
        aa=abs(Erf(end-1,:))./abs(Erf(end-2,:));
        bb=abs(Erf(end-2,:))./abs(Erf(end-3,:));
        krNext2=2.*dkr2(2,:).*aa-dkr2(1,:).*bb;          
        Erf(end,:)=Erf(end-1,:).*exp(1i.*krNext2.*mesh.dr);        
end
%%
% [r_smoothfct,f_smoothfct]=get_smoothfunction(mesh,Erf,'supergaussian');
% Erf=abs(Erf).*f_smoothfct.*exp(1i.*unwrap(angle(Erf)).*f_smoothfct);
if strcmp(globproperties.mode,'debug')
    subplot(3,2,[3 4])
    plot(mesh.f,[(angle(Erf(1:3,:)))]); xlim([-1e14 1e14])
    subplot(3,2,[5 6])
    plot(mesh.f,[(angle(Erf(end-2:end,:)))]); xlim([-1e14 1e14])
    pause(0.1);
end

end