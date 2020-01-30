%Plot stuff
% set(0,'defaultAxesFontSize',10)
%% %%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,[1 2])
plot(1e6.*mesh.r(2:end),real(log(Erf(2:end,pulse.pfmid)./Erf(1:end-1,pulse.pfmid))))
xlabel('radius r (µm)')
ylabel('ln[ E(r_i) / E(r_{i-1}) ]')
set(gca,'YTickLabel',[])
title('real part of dE(r_i,\omega_0) ')
legend('const. Div.','open')


subplot(2,2,[3 4])
plot(1e6.*mesh.r(2:end),imag(log(Erf(2:end,pulse.pfmid)./Erf(1:end-1,pulse.pfmid))))
xlabel('radius r (µm)')
ylabel('ln[ E(r_i) / E(r_{i-1}) ]')
legend('const. Div.','open','LOCATION','East')
title('imaginary part of dE(r_i,\omega_0) ')
set(gca,'YTickLabel',[])



%% %%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(1e6.*mesh.r,[matTprop(:,1),matTprop(:,end-2)])
title('SPM+Plasma')
xlabel('radius r (µm)')
ylabel('|E(r)|^2 (arb.u.)')
set(gca,'YTickLabel',[])
legend('z=0','z=0.3mm (const. Div.)','z=0.3mm (open BC)')
set(0,'defaultAxesFontSize',10)
%% %%%%%%%%%%%%%%%%%%%%%%%%%
figure;
myz=0:mesh.dz:zprop;
surf(1e3.*myz(1:end-30),1e6.*mesh.r,matTprop(:,1:end-30))
ylabel('radius r (µm)')
xlabel('propagation z (mm)')
title('SPM+Plasma, open BC')






%% %%%%%%%%%%%%%%%%%%%%%%%%%
figure;
set(0,'defaultAxesFontSize',10)
imagesc(mesh.f.*1e-12,mesh.r.*1e6,abs(Erf).^2)
xlabel('frequency f (THz)')
ylabel('radius r (µm)')
title('Open BC, z=0.6mm')


%% %%%%%%%%%%%%%%%%%%%%%%%%%
myz=0:mesh.dz:zprop;
figure;
plot(1e6.*mesh.r,matTprop(:,end-50))
bound=myz(end-50);
bound2=find(myz>bound,1);
hold on
plot(1e6.*mesh.r(1:89),matTprop(1:89,bound2))
legend('open','const. Div.')
ylabel('|E(r)|^2 arb.u.')
title('SPM+Plasma, z=0.85mm')
%% Test Filters

Ei=abs(pulse.Ert(1,:)).^2./max(abs(pulse.Ert(1,:)).^2);
GaussDiv=calc_supergaussian(mesh.t,pulse.t_Ie2*1,10,0);
figure; plot(mesh.t,[Ei;mesh.Gfilter_T;GaussDiv]);

Ei=abs(Er(1,:)).^2./max(abs(Er(1,:)).^2);
GaussDiv=calc_supergaussian(mesh.t,pulse.t_Ie2,10,0);
figure; plot(mesh.t,[Ei;mesh.Gfilter_T;GaussDiv]);

Efi=abs(myfft(sumpulse.Ert(1,:),mesh)).^2./max(abs(myfft(sumpulse.Ert(1,:),mesh)).^2);
figure; plot(mesh.f,[Efi;mesh.Tfilter_LR;mesh.Gfilter_T])




