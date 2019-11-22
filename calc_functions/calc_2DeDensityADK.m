%Calculates W_adk for a given E field and a given Gas
function [prod_adk]=calc_2DeDensityADK(prod_adk,mesh,medium,beam,pulse)
% Ert=abs(0.5.*((Ert.*exp(pulse.carrier))+conj(Ert.*exp(pulse.carrier))));
prod_adk=abs(real(prod_adk));
smooth=calc_supergaussian(mesh.t,2*pulse.t_Ie2,10,0);
prod_adk=smooth.*prod_adk;
%Quantum Numbers
l=medium.l;
m=medium.m;
Z=1;                                                                       %First ionization
f=(2*l+1)*factorial(l+abs(m))/(2^abs(m)*factorial(abs(m))*factorial(l-m));
%% Conversion into Atomic Units
rbohr=4*pi*const.eps0*const.hbar^2/(const.m_e*const.e^2);   
E_hartree=const.m_e*const.e^4/(4*const.eps0^2*const.h^2);                  %Hartree Energy [J]
Efield_hartree=E_hartree/(rbohr*const.e);                                  %Efield from Hartree Energy [V/m]
Eg_au=medium.Eg/E_hartree;                                                        %Energy into atomic units a.u.
prod_adk=prod_adk./Efield_hartree;                                                   %Field into atomic untis a.u.  
dt_atu=mesh.dt./(const.hbar/E_hartree);                                    %Time into atomic time units a.t.u.

%% Calc ionization density n_e with ADK model
n_star=Z*(1/sqrt(2*Eg_au));                                                %Effective principal quantum number
E_0_au=(2*Eg_au)^(3/2);

% fraction=E_0_au./(Ert);
% fraction(isnan(fraction))=0;
% 
% w_factor=(Ert).^(1.5-2*n_star);
% w_factor(isinf(w_factor))=0;

% prod=w_factor.*exp(-(2/3).*fraction);
prod_adk=((prod_adk).^(1.5-2*n_star)).*exp(-(2/3).*E_0_au./(prod_adk));
prod_adk(isinf(prod_adk))=0;
prod_adk(isnan(prod_adk))=0;

prod_adk=f/(8*pi*n_star).*((4*exp(1)*E_0_au./(n_star)).^(2*n_star)).*sqrt(3./(pi*E_0_au*(2*Eg_au))).*prod_adk;  %Ionization rate
% P_adk=1-exp(-dt_atu.*W_adk);                                             %Ioniz. propability
prod_adk=medium.n_gas.*(1-exp(cumsum(-dt_atu.*prod_adk,2)));                      %Electron density



end