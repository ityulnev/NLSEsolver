%Calculate propagated electric field after distance Lz via fourier split
%step
%d/dz(E) from SVEA in 1D
%Implemented Effects: GVD, SPM, Attenuation, Ionization
function [Erf,matTprop]=do_adaptiveFourierSplitStep2D(mesh,beam,medium,pulse,boundary)

rayl=calc_zrayleigh(beam,mesh,pulse,0);
wprop=calc_wprop(rayl.win,rayl.zr,rayl.zr);

%Build vector-analogon to finite difference matrix M_fd
M_fd=M_FiniteDifference_init(mesh,boundary);

Erf=pulse.Erf;
matTprop=trapz(mesh.t,medium.Iconst.*abs(pulse.Ert).^2,2);
% counter=0;
zprop=0;
dh=mesh.dz;


while zprop<rayl.zr*3
myerrorZ=1;  
    while myerrorZ>0.01
        %SST+SPM+DIV+ION via Runge Kutta
        [Erf1]=do_LowStoreRK(mesh,pulse,beam,medium,Erf,M_fd,dh);
        Erf2=Erf;
        for mm=1:2
            [Erf2]=do_LowStoreRK(mesh,pulse,beam,medium,Erf2,M_fd,dh/2);
        end
        %Error estimated from maximum value
        max1=max(max(abs(Erf1),[],1),[],2);
        max2=max(max(abs(Erf2),[],1),[],2);
%         energy1=trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,abs(Erf1).^2,2),1);
%         energy2=trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,abs(Erf2).^2,2),1);
        myerrorZ=abs(max1-max2)/max1;
        epsilon=1e-6;
        if myerrorZ<epsilon*0.01
            Erf=Erf1;
            zprop=zprop+dh;
            dh=dh*2;
        elseif myerrorZ>epsilon
            dh=dh/2;
        else
            Erf=Erf1;  
            zprop=zprop+dh;
        end
        if dh<1e-9
           dh=1e-9; 
        dd=errordlg('Error in my stepsize','Warning');
        uiwait(dd)
        end
        
    end
%             counter=counter+1;
%             if counter==100
%             matTprop=[matTprop,trapz(mesh.t,medium.Iconst.*abs(Ert).^2,2)];
%             counter=0;
% %             plot(mesh.r,[matTprop(:,1),matTprop(:,end)]); legend('Initial','Propagated')
%             Q=2*pi.*trapz(mesh.r,transpose(mesh.r).*matTprop(:,end));
%             condition=abs(beam.Q_In-2*Q)/abs(beam.Q_In)>0.1;
%             test_errorMSG(condition,'Error in Pulse energy')
%             end
            %% plot
            Qout=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,medium.Iconst.*abs(Erf).^2,2),1);
%             Ert=myifft(Erf,mesh);
%             plot(mesh.r,[abs(pulse.Ert(:,mesh.tmid)).^2,abs(Ert(:,mesh.tmid)).^2])
%             plot(mesh.f,[abs(pulse.Erf(end,:)).^2;abs(Erf(end,:)).^2]); xlim([0e14 8e14]) 
%             plot(mesh.t,[abs(pulse.Ert(1,:)).^2;abs(Ert(1,:)).^2]); xlim([-1e-13 1e-13])
%             imagesc(mesh.f,mesh.r,abs(Erf).^2); xlim([1e14 6e14])       
            imagesc(mesh.t,mesh.r,abs(myifft(Erf,mesh)).^2); xlim([-1e-13 1e-13])   
            title(['dh=',num2str(dh*1000),'mm','  ','z=',num2str(zprop*1000),'mm','  ','Qout=',num2str(Qout*1000),'mJ']);
            pause(0.1);             
%             if strcmp(lastwarn,'myErrorID')                
%                 save('mysavefromError2.mat');
%                 quit;
%             end
end

end