% Calculates part of Equation in time domain in optimized form to save
% memory
%Eout(t)=f(Ein)=[SPMfocusing+PLSMdefocusing]*Et+[IONloss]*Et
% SPM= const*abs(E)^2
% PLSM= const*n_e(E) 
% ION= const*d/dt[n_e(E)] *1/Re[E]^2
function [Ert]=calc_mainFctOptimizeTime(beam,mesh,medium,pulse,Ert)
%% Self Phase Modulation
const_SPM=-(medium.n2/const.c)*medium.Iconst;
%% Ionization with ADK model // Energy loss via Ionization
const_ION=-medium.Eg/(2*const.eps0*const.c*medium.n0);
%% Plasma Defocusing
% const_PLSM=(1i*pulse.w0/const.c)*const.e^2/((const.eps0*const.m_e)*(2*pulse.w0.^2));
const_PLSM=-const.e^2/(2*const.c*medium.n0*const.m_e*const.eps0);
%% Optimized for less memory use
Ert=handle_NaNInf(gradient(const_SPM.*abs(Ert).^2.*Ert,mesh.dt))+const_PLSM.*cumsum(calc_2DeDensityADK(Ert,mesh,medium,beam,pulse).*mesh.dt.*Ert,2)+handle_NaNInf(const_ION.*(gradient(calc_2DeDensityADK(Ert,mesh,medium,beam,pulse),mesh.dt)./(real(Ert).^2)).*Ert);
Ert(isnan(Ert))=0;%0/0 is NaN ** @zero intensity all should be zero!
Ert(isinf(Ert))=0;
end

