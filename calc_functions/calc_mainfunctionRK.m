% Calculate f'=f(z,E) for Runge Kutta
function [Erf]=calc_mainfunctionRK(mesh,pulse,beam,medium,Erf,M_fd)
%%
Et=myifft(Erf,mesh);
%% Self Phase Modulation // Self Steepening
SPM=-1i*(2*pi.*beam.f0.*medium.n2/const.c).*medium.Iconst.*abs(Et).^2;
% SPMSST=(2*pi.*(mesh.f))./(pulse.w0).*myfft((SPM.*Et),mesh);%
% SPMSST(abs(SPMSST)<1e-10.*max(max(abs(SPMSST))))=0;%take care of numerical errors which might occur            
%% Ionization with ADK model // Energy loss via Ionization
[n_e,Eg,n_gas]=calc_2DeDensityADK(Et,mesh,medium,beam,pulse);
dNedt=diff(n_e,1,2)./mesh.dt;
dNedt=[dNedt,zeros(mesh.rlength,1)];
ION=-Eg.*dNedt./(2.*medium.Iconst.*abs(Et).^2);
ION(isnan(ION))=0;%0/0 is NaN ** @zero intensity all should be zero!
%% Plasma Defocusing
wp2=const.e^2.*n_e./(const.eps0*const.m_e);
PLSM=(1i*pulse.w0/const.c).*wp2./(2.*pulse.w0.^2);
NL=myfft((SPM+PLSM+ION).*Et,mesh);
%% Group velocity dispersion
% GVD=-1i.*(2*pi.*mesh.f-pulse.w0).^2.*medium.kGVD(pulse.pfmid)./2;
%% Divergence in cylinder coordinates
[Etransv]=do_2Dfinitedifference(mesh,medium,Erf,M_fd);
% fboundslr=find_mybounds(mesh.f,-beam.f0/2,beam.f0/2);
% cutk=[zeros(1,fboundslr(1,1)-1),medium.k(fboundslr(1,1):fboundslr(1,2)),zeros(1,mesh.flength-(fboundslr(1,2)))];
DIV=(-1i./(2.*(medium.k0))).*Etransv;
DIV(isinf(DIV))=0;
DIV(isnan(DIV))=0;
%%
Erf=NL+DIV;
end
