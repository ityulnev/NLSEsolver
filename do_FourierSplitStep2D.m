%Calculate propagated electric field after distance Lz via fourier split
%step
%d/dz(E) from SVEA in 1D
%Implemented Effects: GVD, SPM, Attenuation, Ionization
function [Erf,matTprop]=do_FourierSplitStep2D(mesh,beam,medium,pulse,boundary)

rayl=calc_zrayleigh(beam,mesh,pulse,0);
wprop=calc_wprop(rayl.win,rayl.zr,rayl.zr);

%Build vector-analogon to finite difference matrix M_fd
M_fd=M_FiniteDifference_init(mesh,boundary);

Erf=pulse.Erf;
matTprop=trapz(mesh.t,medium.Iconst.*abs(pulse.Ert).^2,2);
counter=0;
    for m=1:(mesh.zlength)

            %SST+SPM+DIV+ION via Runge Kutta
            [Erf]=do_LowStoreRK(mesh,pulse,beam,medium,Erf,M_fd,m);
            Erf(isnan(Erf))=0;
            test_ArrayValues(Erf);
            Ert=myifft(Erf,mesh);        
            
            counter=counter+1;
            if counter==100
            matTprop=[matTprop,trapz(mesh.t,medium.Iconst.*abs(Ert).^2,2)];
            counter=0;
%             plot(mesh.r,[matTprop(:,1),matTprop(:,end)]); legend('Initial','Propagated')
            end
            %% plot
%             Qout=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,fiber.Iconst.*abs(Erf).^2,2),1);
%             plot(mesh.r,[abs(pulse.Ert(:,mesh.tmid)).^2,abs(Ert(:,mesh.tmid)).^2])
%             plot(mesh.f,[abs(pulse.Erf(1,:)).^2;abs(Erf(1,:)).^2]); xlim([3e14 5e14]) 
%             plot(mesh.t,[abs(pulse.Ert(1,:)).^2;abs(Ert(1,:)).^2]); xlim([-1e-13 1e-13])
%             imagesc(mesh.f,mesh.r,abs(Erf).^2); xlim([1e14 6e14])       
%             imagesc(mesh.t,mesh.r,abs(myifft(Erf,mesh)).^2); xlim([-1e-13 1e-13])   
%             title(['z=',num2str(m*mesh.dz),'Qout=',num2str(Qout)]);
%             pause(0.1);
% 
    end
end

