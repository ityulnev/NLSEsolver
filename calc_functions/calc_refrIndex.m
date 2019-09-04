%Calculates the refractive index of a certain gas for a specific wavelength
%according to https://refractiveindex.info/?shelf=main&book=Ne&page=Bideau-Mehu
function [n,ref_n]=calc_refrIndex(wavelength,gas,pressure,T)

wavelength=wavelength.*1e6; %####FORMULA TAKES INPUT IN [um]####
switch gas
    case 'Neon'
    %Refractive Index at T=273.14[K] and 1.01325[bar] (101325Pa)
    ref_n=1+0.00128145./(184.661-wavelength.^(-2))+0.0220486./(376.84-wavelength.^(-2)); %https://refractiveindex.info/?shelf=main&book=Ne&page=Bideau-Mehu
    case 'Neon_n2'
    %n0+n2*Ipeak for pressure scaling of n2     
    ref_n=1+0.00128145./(184.661-wavelength.^(-2))+0.0220486./(376.84-wavelength.^(-2));
    ref_n=ref_n+7.5e-25*2.224e18;
end

%calculate Polarizability alpha for reference conditions with https://de.wikipedia.org/wiki/Clausius-Mossotti-Gleichung
%alpha: property of specific molecule and should be a constant in ideal Gas
Rgas=const.kb*const.Na;
ref_T=273.15;%[K]
ref_pressure=1.01325;
alpha=3*const.eps0*Rgas*ref_T./((3./(ref_n.^2-1)+1).*const.Na*ref_pressure);

%calculate pressure dependent refractive Index of the gas
X=(3*const.eps0*Rgas*T)./(alpha.*const.Na*pressure);
eps_r=1+3./(X-1);%permittivity of gas
n=sqrt(eps_r);

end