%Sets boundary conditions left and right of radial coordinate E(r)
%default = "open" is linear extrapolation to r(1) and r(end)
function [Er]=set_boundaries(mesh,pulse,medium,Er,M_fd,mm,Eprev)
%% Left Boundary
switch M_fd.boundary{1,1}
    case 'constDiv'
        Er(1,:)=Er(1,:);
    case 'const'
        Er(1,:)=Er(2,:);
    case 'open'
        [~,fleft_smoothfct]=get_smoothfunction(mesh,pulse,Er(2,:),'supergaussian');
%         fleft_smoothfct=calc_supergaussian(mesh.t,pulse.t_Ie2,10,0);
        dE=handle_NaNInf(log(Er(3:4,:)./Er(2:3,:))./(1i*mesh.dr));
        kr_LNext=2.*dE(1,:)-1.*dE(2,:);
        kr_LNext=kr_LNext.*fleft_smoothfct;
        Er(1,:)=Er(2,:).*exp(-1i.*kr_LNext.*mesh.dr);
 
end
%% Right Boundary
switch M_fd.boundary{1,2}
    case 'const'
        Er(end,:)=Er(end,:);
    case 'open'
        dkr=handle_NaNInf(log(Er(end-2:end-1,:)./Er(end-3:end-2,:))./(1i*mesh.dr));
        krNext=2.*dkr(2,:)-dkr(1,:);                                          %delKr(Nr+1)=delKr(Nr)+dKr, with dKr=delKr(Nr)-delKr(Nr-1), dKr is a constant value! b.c. linear function with dr->const
        [~,fright_smoothfct]=get_smoothfunction(mesh,pulse,Er(end,:),'supergaussian');
%         fright_smoothfct=calc_supergaussian(mesh.t,pulse.t_Ie2,10,0);
        krNext=krNext.*fright_smoothfct;
        Er(end,:)=Er(end-1,:).*exp(1i.*(krNext).*mesh.dr);                             
end
%% Plot Boundary Conditions
if strcmp(globproperties.mode,'debug')
    subplot(4,2,[3 4])
    plot(mesh.t,[(angle(Er(1:3,:)))]); xlim([-100e-15 100e-15])
    subplot(4,2,[5 6])
    plot(mesh.t,[(angle(Er(end-2:end,:)))]); xlim([-100e-15 100e-15])
    pause(0.1);
end

end