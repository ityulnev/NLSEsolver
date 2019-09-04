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

figure; imagesc(matTprop(:,1:end))
figure; surf(matTprop(:,1:end))
figure; plot(mesh.r',[norm_fields(matTprop(:,1)',matTprop(:,end)','indiv')])
figure; plot(mesh.r,[matTprop(:,1),matTprop(:,end)])


figure; imagesc(abs(Erf).^2)
figure; plot(mesh.f,[abs(pulse.Erf(1,:)).^2;abs(Erf(1,:)).^2])

Ert=myifft(Erf,mesh);
figure; imagesc(abs(Ert).^2)
figure; plot(mesh.f,[abs(pulse.Ert(1,:)).^2;abs(Ert(1,:)).^2;abs(Ert(300,:)).^2])
figure; plot(mesh.f,[(pulse.Ert(1,:));abs(pulse.Ert(1,:));(Ert(1,:));abs(Ert(1,:))])


figure; plot(mesh.f,[norm_fields(abs(pulse.Erf(1,:)),abs(Erf(1,:)),'indiv')])
%%
set(gca, 'YScale', 'log')

E1=gaussfilter.*(E_opt);
E2=(E_opt);
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
