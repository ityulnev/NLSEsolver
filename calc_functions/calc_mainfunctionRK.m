% Calculate f'=f(z,E) for Runge Kutta
function [Erf]=calc_mainfunctionRK(mesh,pulse,beam,medium,Erf,M_fd)
%%
% Et=myifft(Erf,mesh);
% %% SPM+SST
% SPM=-1i*(2*pi.*beam.f0.*medium.n2/const.c).*medium.Iconst.*abs(Et).^2;
% % SPMSST=(2*pi.*(mesh.f))./(pulse.w0).*myfft((SPM.*Et),mesh);%
% % SPMSST(abs(SPMSST)<1e-10.*max(max(abs(SPMSST))))=0;%take care of numerical errors which might occur            
% %% Ionization with ADK model
% [n_e,Eg]=calc_2DeDensityADK(Et,mesh,medium);
% % dNedt=diff(n_e,1,2)./mesh.dt;
% % dNedt=[dNedt,zeros(mesh.rlength,1)];
% % ION=-Eg.*dNedt./(2.*medium.Iconst.*abs(Et).^2);
% % ION(isnan(ION))=0;%0/0 is NaN ** @zero intensity all should be zero!
% %% Plasma Defocusing
% wp2=const.e^2.*n_e./(const.eps0*const.m_e);
% PLSM=(1i*pulse.w0/const.c).*wp2./(2.*pulse.w0.^2);
% NL=myfft((10.*SPM+5.*PLSM).*Et,mesh);
%% Group velocity dispersion
% GVD=-1i.*(2*pi.*mesh.f-pulse.w0).^2.*medium.kGVD(pulse.pfmid)./2;
%% Divergence in zylinder coordinates
[Etrans]=do_2Dfinitedifference(mesh,medium,Erf,M_fd);
DIV=(-1i./(2.*(medium.k+medium.k0))).*Etrans;
% DIV(:,1:find(mesh.f>beam.f0*0.20,1))=0; %1/k(f~0) leads to big errors! 
% DIV(abs(DIV)<max(max(abs(DIV),[],1),[],2).*1e-6)=0;
DIV(isinf(DIV))=0;
DIV(isnan(DIV))=0;
%% Take care of Error with smoothing function
Imax=max(abs(Erf(1,:)).^2);
smoothbound=find(abs(Erf(1,:).^2)>Imax.*1e-9,1);
halfwidth=mesh.df*abs(mesh.fbound-smoothbound);
smoothfct=calc_supergaussian(mesh.f,halfwidth,10);
DIV=DIV.*smoothfct;

%%
NL=0;
Erf=NL+DIV;
end
