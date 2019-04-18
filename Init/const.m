classdef const
    properties (Constant)
    c=3e8;                                                                 %speed of light [m/s]
    eps0=8.854e-12;                                                        %vacumm permitivity [J/(V^2*m)]
    h=6.626e-34;                                                           %Plank Constant [Js]
    m_e=9.1094e-31;                                                        %Electron mass [kg]
    hbar=6.626e-34/(2*pi);                                                 %Plank Constant [Js/rad]
    e=1.602e-19;                                                           %Charge of Electron [C]
%     rbohr=4*pi*eps0*hbar^2/(m_e*e^2);                                      %Bohrradius
    Na=6.02214e23;                                                         %Avogadro Constant[1/mol]
    kb=1.380649e-23;                                                        %Boltzmann Constant[J/K]
    end    
end