% Calculate f'=f(z,E) for Runge Kutta
function [result]=calc_mainfunctionRK(mesh,pulse,beam,medium,Erf,M_fd)
%%
% Et=myifft(Erf,mesh);
% % SPM+SST
% SPM=-1i*(2*pi.*beam.f0.*fiber.n2/const.c).*fiber.Iconst.*abs(Et).^2;
% % SPMSST=(2*pi.*(mesh.f))./(pulse.w0).*myfft((SPM.*Et),mesh);%
% % SPMSST(abs(SPMSST)<1e-10.*max(max(abs(SPMSST))))=0;%take care of numerical errors which might occur            
% %% Ionization with ADK model
% [n_e,Eg]=calc_2DeDensityADK(Et,mesh,fiber);
% % dNedt=diff(n_e,1,2)./mesh.dt;
% % dNedt=[dNedt,zeros(mesh.rlength,1)];
% % ION=-Eg.*dNedt./(2.*fiber.Iconst.*abs(Et).^2);
% % ION(isnan(ION))=0;%0/0 is NaN ** @zero intensity all should be zero!
% %% Plasma Defocusing
% wp2=const.e^2.*n_e./(const.eps0*const.m_e);
% PLSM=(1i*pulse.w0/const.c).*wp2./(2.*pulse.w0.^2);
% %%

% wplasma=sqrt(n_e.*const.e^2./(const.eps0*const.m_e));
% plasma_defocusing=-(1/(2*const.c)).*trapz(mesh.t,wplasma.^2.*Ert,2);
%             
%Group velocity dispersion
% % GVD=-1i.*(2*pi.*mesh.f-pulse.w0).^2.*medium.kGVD(pulse.pfmid)./2;    
%Divergence in zylinder coordinates
[Etrans]=do_2Dfinitedifference(mesh,medium,Erf,M_fd);
Etrans(abs(Etrans).^2<max(max(abs(Etrans).^2)).*1e-9)=0;
% NL=myfft((0.*SPM+0.*PLSM).*Et,mesh);
NL=0;
DIV=(-1i./(2.*medium.k)).*Etrans;
DIV(isinf(DIV))=0;
DIV(isnan(DIV))=0;
% div=0;
result=NL+DIV;

%plot(mesh.r,[unwrap(angle(spm(:,pulse.pfmid)))])
%plot(mesh.r,[unwrap(angle(spm(:,pulse.pfmid))),unwrap(angle(div(:,pulse.pfmid)))])
end
