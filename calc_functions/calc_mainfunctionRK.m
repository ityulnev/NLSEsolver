% Calculate f'=f(z,E) for Runge Kutta
function [Ert]=calc_mainfunctionRK(mesh,pulse,medium,Ert,M_fd)
% Group velocity dispersion
const_GVD=-1i.*(2*pi.*mesh.f-pulse.w0).^2.*medium.k2_w0./2;

% full simulation
Ert=handle_NaNInf(myifft(mesh.Gfilter_T.*const_GVD.*myfft(abs(Ert).*exp(1i.*pulse.w0.*mesh.t.*ones(mesh.rlength,1)),mesh),mesh))+(calc_mainFctOptimizeTime(mesh,medium,pulse,Ert))+(const.c/(2.*medium.n0)).*do_2Dfinitedifference(mesh,medium,cumsum(Ert.*mesh.dt,2),M_fd);
% no transverse effect!
% Ert=handle_NaNInf(myifft(const_GVD.*myfft(abs(Ert).*exp(pulse.carrier),mesh),mesh))+(calc_mainFctOptimizeTime(beam,mesh,medium,pulse,Ert));

Ert=(Ert).*mesh.Gfilter_T.*mesh.Gfilter_R';
Ert=do_filter(Ert,'TanhandGaussian','inF',mesh);
end
