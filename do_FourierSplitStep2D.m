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
h=mesh.dz;
zprop=0;
    while zprop<(rayl.zr*3)
            %SST+SPM+DIV+ION via Runge Kutta
            [Erf]=do_LowStoreRK(mesh,pulse,beam,medium,Erf,M_fd,h);
            condition1=sum(sum(isinf(Erf),1),2)>=1;
            test_errorMSG(condition1,'Error in Erf: Inf values')
            condition2=sum(sum(isnan(Erf),1),2)>=1;
            test_errorMSG(condition2,'Error in Erf: NaN values')

            counter=counter+1;
            zprop=zprop+h;
            if counter==1
                Ert=myifft(Erf,mesh);   
                matTprop=[matTprop,trapz(mesh.t,medium.Iconst.*abs(Ert).^2,2)];
                condition3=sum(sum(isinf(matTprop),1),2)>=1;
                test_errorMSG(condition3,'Error in Ert: Inf values')
                counter=0;
                Q=2*pi.*trapz(mesh.r,transpose(mesh.r).*matTprop(:,end));
                dQ=abs(pulse.Energy-Q)/abs(pulse.Energy);
%                 condition=dQ>0.1;
%                 test_errorMSG(condition,'Error in Pulse energy')
                plot(mesh.r,[matTprop(:,1),matTprop(:,end)]); legend('Initial','Propagated');
                title(['z=',num2str(zprop*1000),'mm','  ','Qout=',num2str(Q*1000),'mJ','  ','dQ=',num2str(dQ)]);
                pause(0.1);
            end                    
            %% plot
%             Qout=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,medium.Iconst.*abs(Erf).^2,2),1);
%             plot(mesh.r,[abs(pulse.Ert(:,mesh.tmid)).^2,abs(Ert(:,mesh.tmid)).^2])
%             plot(mesh.f,[abs(pulse.Erf(end,:)).^2;abs(Erf(end,:)).^2]); xlim([3e14 5e14]) 
%             plot(mesh.t,[abs(pulse.Ert(1,:)).^2;abs(Ert(1,:)).^2]); xlim([-1e-13 1e-13])
%             imagesc(mesh.f,mesh.r,abs(Erf).^2); xlim([1e14 6e14])       
%             imagesc(mesh.t,mesh.r,abs(myifft(Erf,mesh)).^2); xlim([-1e-13 1e-13])   
%             title(['z=',num2str(zprop)]);
%             pause(0.1)
            if strcmp(lastwarn,'myErrorID')
                        dd=errordlg('Stopped due to Error','Warning');
                        uiwait(dd)
%                 save('mysavefromError2.mat');
%                 quit;
            end
% 
    end
end

