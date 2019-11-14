%Calculate propagated electric field after distance Lz via Rk with FT
%d/dz(E) from SVEA in 1D
%Implemented Effects: GVD, SPM, Attenuation, Ionization(Defocusing & Loss)
function [Er,E_tmax,E_t0max,ne_tmax,ne_t0max,zprop,Qhist,whist,IonizLvl]=do_WaveEqSolver(mesh,beam,medium,pulse,Er,boundary,comment)

rayl=calc_zrayleigh(beam,mesh,medium,pulse,0);
wprop=calc_wprop(rayl.win,rayl.zr,rayl.zr);

%Build vector-analogon to finite difference matrix M_fd
M_fd=M_FiniteDifference_init(mesh,boundary);

matTprop=trapz(mesh.t,medium.Iconst.*abs(Er).^2,2);
counter=0;
h=mesh.dz;
zprop=0;
Qhist=pulse.Energy;
m=0;
whist=mesh.dr*2+(find(abs(pulse.Ert(:,pulse.ptmid)).^2<abs(pulse.Ert(1,pulse.ptmid))^2/exp(2),1)-1)*mesh.dr;
[n_e]=calc_2DeDensityADK(Er,mesh,medium,beam,pulse);
IonizLvl=max(n_e(1,:))/medium.n_gas;
n_ein=n_e;

%Initial matrizes
tmax= max(real(pulse.Ert).^2,[],2)==real(Er).^2;
E_tmax=Er(tmax);
ne_tmax=n_e(tmax);
E_t0max=Er(:,mesh.indexfmid);
ne_t0max=n_e(:,mesh.indexfmid);

    while zprop<(mesh.Lz)
%%        
        m=m+1;
            %SST+SPM+DIV+ION via Runge Kutta
            [Er]=do_LowStoreRK(mesh,pulse,beam,medium,Er,M_fd,h);
            condition1=sum(sum(isinf(Er),1),2)>=1;
            test_errorMSG(condition1,'Error in Erf: Inf values')
            condition2=sum(sum(isnan(Er),1),2)>=1;
            test_errorMSG(condition2,'Error in Erf: NaN values')

            counter=counter+1;
            zprop=m*mesh.dz;
            if counter==1
%                 Ert=myifft(Erf,mesh); 
%                 matTprop=[matTprop,trapz(mesh.t,medium.Iconst.*abs(Er).^2,2)];
%                 condition3=sum(sum(isinf(matTprop),1),2)>=1;
%                 test_errorMSG(condition3,'Error in Ert: Inf values')

                % Ionization level
                [n_e]=calc_2DeDensityADK(Er,mesh,medium,beam,pulse);
                IonizLvl=[IonizLvl,max(n_e(1,:))/medium.n_gas];
                % Save Er at t0 and at Real(E).^2 max
                tmax= max(real(Er).^2,[],2)==real(Er).^2;
                
                E_tmax=[E_tmax,Er(tmax)];
                ne_tmax=[ne_tmax,n_e(tmax)];
                E_t0max=[E_t0max,Er(:,mesh.indexfmid)];
                ne_t0max=[ne_t0max,n_e(:,mesh.indexfmid)];
                                
                counter=0;
                Q=2*pi.*trapz(mesh.r,transpose(mesh.r).*matTprop(:,end));
                dQ=(pulse.Energy-Q)/abs(pulse.Energy);
                Qhist=[Qhist,Q];
                waist=mesh.dr*2+(find(abs(Er(:,pulse.ptmid)).^2<abs(Er(1,pulse.ptmid))^2/exp(2),1)-1)*mesh.dr;
                whist=[whist,waist];
                if strcmp(globproperties.mode,'debug')
                    subplot(4,2,[1 2])
%                     plot(mesh.r,[matTprop(:,1),matTprop(:,end)]); legend('Initial','Propagated');
                    plot(mesh.r,[real(E_tmax(:,1)),real(E_tmax(:,end))]); legend('Initial','Propagated');
                    title(['z=',num2str(zprop*1000),'mm','  ','Qout=',num2str(Q*1000),'mJ','  ','dQ=',num2str(dQ)]);
                    subplot(4,2,[7 8])
                    plot(mesh.t,abs(Er(1,:))); xlim([-5e-14 5e-14])
                    pause(0.1);
                end

            end                    
            
            
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

