%%
pp3=zeros(1,30);
pp2=0;
figure;
for hh=35900:5:35930

pp=polyfit(mesh.r(2:30)',angle(Erf(2:30,hh)),3);
pp2=polyval(pp,mesh.r(1:30));
pp3=[pp3;pp2];
hold on
plot(mesh.r(1:30),pp2)
pause(0.1)
end
figure;
plot(mesh.f(35900:5:35930),[pp3(2:end,1)])

hold on; plot(mesh.f(35900:5:35930),angle(Erf(1:4,35900:5:35930)))
%% Debug plots

figure; imagesc(matTprop(:,1:end-3))
figure; surf(matTprop(:,1:end-3))
figure; plot(mesh.r',[norm_fields(matTprop(:,1)',matTprop(:,end)','indiv')])
figure; plot(mesh.r,[matTprop(:,1),matTprop(:,end-3)])


figure; imagesc(abs(Erf).^2)
figure; plot(mesh.t,[norm_fields(abs(pulse.Ert(1,:)).^2,abs(Erf(1,:)).^2,'indiv')])

% Ert=myfft(Erf,mesh);
Ert=Er;
figure; imagesc(abs(Ert).^2)
figure; plot(mesh.f,[abs(pulse.Ert(1,:)).^2;abs(Ert(1,:)).^2;abs(Ert(end,:)).^2])
figure; plot(mesh.f,[(pulse.Ert(1,:));abs(pulse.Ert(1,:));(Ert(1,:));abs(Ert(1,:))])


figure; plot(mesh.f,[norm_fields(abs(pulse.Erf(1,:)),abs(Erf(1,:)),'indiv')])
%%
Ef=myfft(Er,mesh);
figure; plot(mesh.f,[abs(pulse.Erf(1,:));abs(Ef(1,:))])
figure; imagesc(abs(Ef).^2)
%%
set(gca, 'YScale', 'log')

E1=(calc_mainFctOptimizeTime(beam,mesh,medium,pulse,E_opt));
E2=(gaussfilter.*E1);
Ef2=myfft(E2,mesh);
shift=(mesh.indexfmid-1)/2;
E2sh=[zeros(1,shift),E2(1,1:end-shift)];
Ef2sh=myfft(E2sh,mesh);
Efmirr=fliplr(Ef2);
figure; plot(mesh.f,[abs(E2(1,:));abs(E2m)])
figure; plot(mesh.f,[abs(Ef2(1,:));abs(Ef2sh)])
figure; plot(mesh.f,[unwrap(angle((Ef2(1,:))));unwrap(angle(fftshift(Ef2(1,:))));unwrap(angle(ifftshift(Ef2(1,:))))])

gauss1=calc_supergaussian(mesh.t,pulse.fwhmT.*2,10,0);
gauss2=calc_supergaussian(mesh.f,(pulse.pfmid-LRbounds(1,1)).*mesh.df,10,beam.f0);
E3=E_opt.*gauss1;
E4=myfft(E3,mesh);

Emirr=flip((E2),2);
Eclean=E2-Emirr;
Eclean(:,1:mesh.indexfmid)=0;
figure; plot(mesh.t,[abs(E1(1,:));abs(E3(1,:))])
figure; plot(mesh.f,[filgaussPLSM;abs(E2(1,:))./max(abs(E2(1,:)));abs(Erf(1,:))./max(abs(Erf(1,:)))])
figure; plot(mesh.f,abs(E1(1,:)))

figure; plot(mesh.t,[abs(E1(1,:))./max(abs(E1(1,:)));gauss1])
%%

n_e=calc_2DeDensityADK(Et,mesh,medium,beam,pulse);
figure; imagesc(n_e)
figure; plot(mesh.t,n_e(1,:)./medium.n_gas)
%%

Ett=-1i*(medium.n2/const.c).*medium.Iconst.*abs(E_opt).^2.*E_opt.*gaussfilter;
Eff=(2*pi.*mesh.f).*myfft(Ett,mesh).*filpos;
figure; plot(mesh.t,[abs(Ett(1,:))./max(abs(Ett(1,:)));gaussfilter])
figure; plot(mesh.f,abs(Eff(1,:)).^2)
%%

figure;plot(mesh.t,[norm_fields(unwrap(angle(Erf(end,:))),unwrap(angle(Erf(end-1,:))),krNext,fright_smoothfct,'indiv')])
hold on; plot(mesh.t,norm_fields(unwrap(angle(Erf(end,:))),'indiv'))
hold on; plot(mesh.t,norm_fields(abs(krNext),'indiv'))
figure;plot(mesh.t,[norm_fields(abs(Erf(end,:)),abs(Erf(end-1,:)),krNext,fright_smoothfct,'indiv')])
%% Energy conservation
Nein=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t.*const.c,n_ein,2),1);
Ne=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t.*const.c,n_e,2),1);
Eloss=(Nein-Ne)*medium.Eg
Qhist(1,1)-Qhist(1,2)
Na=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t.*const.c,medium.n_gas.*ones(1,mesh.flength),2),1);

Nemax=2*pi.*trapz(mesh.r,transpose(mesh.r).*n_e(:,end),1);
dNe=Nemax*(mesh.dz*medium.n0pressure/const.c)

% Namax=2*pi.*trapz(mesh.r,transpose(mesh.r).*medium.n_gas.*ones(mesh.rlength,1),1);
Volume=pi*mesh.R^2*(mesh.t(end)*2)*const.c
Ne/Volume
myNe=Na*max(n_e(1,:))./medium.n_gas
Eloss=Nemax.*medium.Eg

Ener=max(n_e(1,:)).*medium.Eg;

const_ION=medium.Eg/(1);%-2*const.eps0*const.c*medium.n0
const_PLSM=const.e^2/(const.m_e);%-2*const.c*medium.n0*const.m_e*const.eps0
% J=const_PLSM.*cumsum(calc_2DeDensityADK(pulse.Ert,mesh,medium,beam,pulse).*mesh.dt.*pulse.Ert,2)+handle_NaNInf(const_ION.*(gradient(calc_2DeDensityADK(pulse.Ert,mesh,medium,beam,pulse),mesh.dt)./(real(pulse.Ert).^2)).*pulse.Ert);
Jin=const_PLSM.*cumsum(calc_2DeDensityADK(pulse.Ert,mesh,medium,beam,pulse).*mesh.dt.*real(pulse.Ert),2);
Jprop=const_PLSM.*cumsum(calc_2DeDensityADK(Er,mesh,medium,beam,pulse).*mesh.dt.*real(Er),2);
RhoEin=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t,(pulse.Ert).*Jin,2),1);
RhoEprop=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t,(Er).*Jprop,2),1);
dE=(abs(RhoEin)-abs(RhoEprop))./(trapz(mesh.r,mesh.r'.*2*pi))



%% Drift Loss @averaged r
velocin=const.e/const.m_e.*cumsum(abs(2*pi.*trapz(mesh.r,mesh.r'.*real(pulse.Ert),1)./trapz(mesh.r,2*pi.*mesh.r')).*mesh.dt);
veloc2=const.e/const.m_e.*cumsum(abs(2*pi.*trapz(mesh.r,mesh.r'.*real(Er),1)./trapz(mesh.r,2*pi.*mesh.r')).*mesh.dt);

Ekin1=const.m_e.*max(velocin).^2./2
Ekin2=const.m_e.*max(veloc2).^2./2

veloc3=const.e/const.m_e.*cumsum(abs(2*pi.*trapz(mesh.r,mesh.r'.*real(Er),1)./trapz(mesh.r,2*pi.*mesh.r')).*mesh.dt);
Ekin3=const.m_e.*max(veloc3).^2./2

(Ekin1-Ekin2)/Ekin1
figure; plot(mesh.t,[velocin;veloc2])
%% Drift Loss @r=1
ind=find(mesh.t>1200e-15,1)-1;
ind=find(max(abs(pulse.Ert(1,:)))./exp(2)<abs(pulse.Ert(1,:)),1,'last')+1
filtr=[ones(1,inde2-1),zeros(1,mesh.flength-(inde2-1))];
velocin=const.e/const.m_e.*cumsum(abs(real(pulse.Ert(1,:))).*mesh.dt);
veloc2=const.e/const.m_e.*cumsum(abs(real(Er(1,:))).*mesh.dt);

Ekin1=const.m_e.*(velocin(ind)).^2./2
Ekin2=const.m_e.*(veloc2(ind)).^2./2
(Ekin1-Ekin2)/Ekin1

mesh.t(ind)
figure; plot(mesh.t,[velocin;veloc2])


EE=trapz(mesh.r,mesh.r'.*trapz(mesh.t,medium.Iconst.*abs(pulse.Ert).^2,2),1).*2*pi
%%
v2=const.e/const.m_e.*cumsum((real(Er)).*mesh.dt);
gamma2=1./sqrt(1-(v2./const.c).^2);
Ekin2=(gamma2-1).*const.m_e.*const.c^2;
Efull2=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t,Ekin2.*n_e,2),1);

dE=(Efull-Efull2)./Efull
(Qhist(1,1)-Qhist(1,end))/Qhist(1,1)

figure; plot(mesh.t,Ekin1./const.e)
figure; plot(mesh.t,(gammaIn-1))
%% Energy loss drift

vIn=(const.e/const.m_e)*2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t,abs(real(pulse.Ert)),2),1)./(2*pi.*trapz(mesh.r,mesh.r',1))
% vIn=(const.e/const.m_e).*trapz(mesh.t,abs(real(pulse.Ert(1,:))),2)

EkinIn=const.m_e/2.*(vIn)^2


%%
[myn_e]=calc_2DeDensityADK(dd.pulse.Ert,dd.mesh,dd.medium,dd.beam,dd.pulse);
IonizLvl=max(myn_e(1,:))/dd.medium.n_gas
%% Self focusing length

Ldf=pi*medium.n0pressure*beam.r_mode^2/beam.wavelength
Lsf=0.367*medium.k0*beam.r_mode^2/sqrt((sqrt(pulse.PpeakTheo/medium.P_crit)-0.852)^2-0.0219)
%% Frequency domain when approaching single cycles
figure; plot(mesh.f.*1e-12,norm_fields([abs(myfft(real(pulse.Ert(1,:)),mesh)).^2],1),'LineWidth',1)
hold on; plot(mesh.f.*1e-12,norm_fields(abs(myfft((pulse.Ert(1,:)),mesh)).^2,1),'LineWidth',1)
hold on; plot(mesh.f.*1e-12,norm_fields(abs(myfft(real(pulse2.Ert(1,:)),mesh)).^2,1),'LineWidth',1)
hold on; plot(mesh.f.*1e-12,norm_fields(abs(myfft((pulse2.Ert(1,:)),mesh)).^2,1),'LineWidth',1)
xlim([-400 400])
xlabel('frequency (THz)')
ylabel('$|\mathcal{F}\{E\}|^2$ (norm.)')
legend('1 cycle','analytic','2 cycle')
my_figure_settings('Lim_FT_Ereal',1)



