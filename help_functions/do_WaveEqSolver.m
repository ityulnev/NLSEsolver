%Calculate propagated electric field after distance Lz via Rk with FT
%d/dz(E) from SVEA in 1D
%Implemented Effects: GVD, SPM, Attenuation, Ionization(Defocusing & Loss)
function [Er,Etrz,zprop,IonizLvl,Zsteps,mm]=do_WaveEqSolver(mesh,beam,medium,pulse,Er,boundary,comment)

rayl=calc_zrayleigh(beam,mesh,medium,pulse,0);
wprop=calc_wprop(rayl.win,rayl.zr,rayl.zr);
Zsteps=10;
%Build vector-analogon to finite difference matrix M_fd
M_fd=M_FiniteDifference_init(mesh,boundary);

matTprop=trapz(mesh.t,medium.Iconst.*abs(Er).^2,2);
counter=1;
h=mesh.dz;
zprop=0;
m=0;
mm=1;
[n_e]=calc_2DeDensityADK(Er,mesh,medium,beam,pulse);
IonizLvl=max(n_e(1,:))/medium.n_gas;
n_ein=n_e;

%shrinked Ert
index_tR=(pulse.indt_Ie2-pulse.ptmid)*2+pulse.ptmid;
index_tL=pulse.ptmid-(pulse.indt_Ie2-pulse.ptmid)*2;
Etrz=Er(:,index_tL:index_tR);


    while zprop<(mesh.Lz)
%%        
        m=m+1;
            %SST+SPM+DIV+ION via Runge Kutta
            [Er]=do_LowStoreRK(mesh,pulse,beam,medium,Er,M_fd,h);
            condition1=sum(sum(isinf(Er),1),2)>=1;
            test_errorMSG(condition1,'Error in Erf: Inf values')
            condition2=sum(sum(isnan(Er),1),2)>=1;
            test_errorMSG(condition2,'Error in Erf: NaN values')

            zprop=m*mesh.dz;
            if counter==Zsteps
%                 Ert=myifft(Erf,mesh); 
%                 matTprop=[matTprop,trapz(mesh.t,medium.Iconst.*abs(Er).^2,2)];
%                 condition3=sum(sum(isinf(matTprop),1),2)>=1;
%                 test_errorMSG(condition3,'Error in Ert: Inf values')

                % Save Er at t0 and at Real(E).^2 max
%                 E_t0max=[E_t0max,Er(:,mesh.indexfmid)];
%                 ne_t0max=[ne_t0max,n_e(:,mesh.indexfmid)];
% %                 Ihist(:,:,1)=Iprev;                
% %                 Ihist(:,:,2)=medium.Iconst.*abs(Er).^2; 
% %                 
% %                 
% %                 Ihist(:,:,3)=medium.Iconst.*abs(Er).^2; 
% %                 alphas=7.8125e-6*(2e-6)^2;                
% %                 dd=diff(Ihist,mesh.dz,3)
% %                 dI=(dd(:,:,1)+dd(:,:,2))./2;
% %                 
% %                 % dk = q*k1 - kq
% %                 zr1=pi*beam.r_mode^2/beam.wavelength;
% %                 k1=(1/zr1)./(1+(myz(m-1)./zr1).^2);
% %                 zrq=q.*pi*beam.r_mode^2/beam.wavelength;
% %                 kq=(1/zrq)./(1+(myz(m-1)./zrq).^2);
% %                 dk.full = q.*(k1)-kq+ alphas.*dI
                    mm=mm+1;
                Etrz(:,:,mm)=Er(:,index_tL:index_tR);

                counter=0;
                if strcmp(globproperties.mode,'debug')
                Q=2*pi.*trapz(mesh.r,transpose(mesh.r).*matTprop(:,end));
                dQ=(pulse.Energy-Q)/abs(pulse.Energy);
                    subplot(4,2,[1 2])
%                     plot(mesh.r,[matTprop(:,1),matTprop(:,end)]); legend('Initial','Propagated');
                    plot(mesh.r,[real(Etrz(:,pulse.ptmid-index_tL,1)),real(Etrz(:,pulse.ptmid-index_tL,mm))]); legend('Initial','Propagated');
                    title(['z=',num2str(zprop*1000),'mm','  ','Qout=',num2str(Q*1000),'mJ','  ','dQ=',num2str(dQ)]);
                    subplot(4,2,[7 8])
                    plot(mesh.t,abs(Er(1,:))); xlim([-5e-14 5e-14])
                    pause(0.1);
                end

            end
            
            counter=counter+1;

            
%% Tests
            if strcmp(lastwarn,'myErrorID')
                switch globproperties.mode
                    case 'debug'
                        dd=errordlg('Stopped due to Error','Warning');
                        uiwait(dd)
                    case 'cluster'
                        %save([date,'savefromerror.mat'],'mesh','pulse','beam','matTprop','Erf','zprop','dQhist','rayl','whist','boundary');
                        save([date,'savefromerror.mat'],'mesh','beam','m','zprop','matTprop','pulse','boundary')
                        quit;
                end
            end
    end
end

