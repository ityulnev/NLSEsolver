% Calculate f'=f(z,E) for Runge Kutta
function [Ert]=calc_mainfunctionRK(mesh,pulse,beam,medium,Ert,M_fd)
%%
% Erf=myifft(Erf,mesh);
gaussfilter=calc_supergaussian(mesh.t,800e-15,10,0);  %mesh.dt.*round(mesh.flength/20)  ||  pulse.t_pulse.*10
%% Self Phase Modulation // Self Steepening
% NL=-1i*(2*pi.*beam.f0.*medium.n2/const.c).*medium.Iconst.*abs(Et).^2;
% SPMSST=(2*pi.*(mesh.f))./(pulse.w0).*myfft((SPM.*Et),mesh);%
% SPMSST(abs(SPMSST)<1e-10.*max(max(abs(SPMSST))))=0;%take care of numerical errors which might occur            
%% Ionization with ADK model // Energy loss via Ionization
% [PLSM]=calc_2DeDensityADK(E_opt,mesh,medium,beam,pulse);
% ION=gradient(ION,mesh.dt);
% ION=-medium.Eg.*ION./(2*const.eps0*const.c*medium.n0.*real(E_opt).^2).*E_opt;
% ION(isnan(ION))=0;%0/0 is NaN ** @zero intensity all should be zero!
% ION(isinf(ION))=0;
%% Plasma Defocusing
% PLSM=cumsum(-const.e^2.*PLSM./(2*const.c*medium.n0*const.eps0*const.m_e).*E_opt.*mesh.dt,2);
% PLSM=(1i*pulse.w0/const.c).*PLSM./(2.*pulse.w0.^2);
% E_opt=-const.e^2/(2*const.c*medium.n0*const.m_e*const.eps0).*cumsum(calc_2DeDensityADK(E_opt,mesh,medium,beam,pulse).*mesh.dt.*E_opt,2);
%% Group velocity dispersion
const_GVD=-1i.*(2*pi.*mesh.f-pulse.w0).^2.*medium.k2_w0./2;
%% Divergence in cylinder coordinates
% filpos=[zeros(1,mesh.indexfmid),ones(1,mesh.flength-mesh.indexfmid)];
% cutk=medium.k_fit;
% cutk(mesh.indexfmid+1:mesh.indexfmid+200)=medium.k_fit(mesh.indexfmid+200);
% Erf=handle_NaNInf((-1i./(2.*(medium.k0))).*do_2Dfinitedifference(mesh,medium,Erf,M_fd))+filpos.*(2*pi.*mesh.f).*myfft(-1i*(medium.n2/const.c).*medium.Iconst.*abs(E_opt).^2.*E_opt.*gaussfilter,mesh);
% fboundslr=find_mybounds(mesh.f,-beam.f0/2,beam.f0/2);
% cutk=[zeros(1,fboundslr(1,1)-1),medium.k(fboundslr(1,1):fboundslr(1,2)),zeros(1,mesh.flength-(fboundslr(1,2)))];
% cutk=medium.k;
% cutk(1:mesh.indexfmid+round((pulse.pfmid-mesh.indexfmid)*(0.20)))=0;
% Erf=Erf;
% Erf(isinf(Erf))=0;
% Erf(isnan(Erf))=0;
%%
% DIV=(const.c/(2.*medium.n0)).*cumsum(do_2Dfinitedifference(mesh,medium,E_opt,M_fd).*mesh.dt,2);
% DIV2=(const.c/(2.*medium.n0)).*do_2Dfinitedifference(mesh,medium,cumsum(E_opt.*mesh.dt,2),M_fd);
%% Time domain effects
% E_opt=myfft(gaussfilter.*PLSM,mesh);
% E_opt=myfft(gaussfilter.*(calc_mainFctOptimizeTime(beam,mesh,medium,pulse,Erf)),mesh);

% Ert=handle_NaNInf(myifft(const_GVD.*myfft(abs(Ert).*exp(pulse.carrier),mesh),mesh))+(calc_mainFctOptimizeTime(beam,mesh,medium,pulse,Ert))+(const.c/(2.*medium.n0)).*do_2Dfinitedifference(mesh,medium,cumsum(Ert.*mesh.dt,2),M_fd);
Ert=handle_NaNInf(myifft(const_GVD.*myfft(abs(Ert).*exp(pulse.carrier),mesh),mesh))+(calc_mainFctOptimizeTime(beam,mesh,medium,pulse,Ert));


% Erf=(const.c/(2.*medium.n0)).*do_2Dfinitedifference(mesh,medium,cumsum(Erf.*mesh.dt,2),M_fd);%(calc_mainFctOptimizeTime(beam,mesh,medium,pulse,E_opt))
% LRbounds=find_bounds(Erf(1,:));
% filgaussPLSM=calc_supergaussian(mesh.f,mesh.df*(LRbounds(1,3)-LRbounds(1,1)),10,mesh.df*(LRbounds(1,2)-mesh.indexfmid));
% filgaussL=[filgaussPLSM(1:LRbounds(1,2)),ones(1,mesh.flength-LRbounds(1,2))];
Ert=(Ert).*gaussfilter;
Ert=do_filter(Ert,'tanhfilterLR','inF',mesh);
%%
% Erf=(E_opt.*filgaussL+Erf);
% Erf=Erf.*mesh.tanhfilterLR;
end
