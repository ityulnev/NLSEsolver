%Create structure with all needed parameters specific to the beam/pulse
classdef beam_init

    properties
    wavelength,f0,w0,Q_In,Q_Out,t_pulse,t_fwhm,r_mode,area_mode,Fluence,alpha,n_cycles
    end
    
    methods 
        function s=beam_init
        s.wavelength=2000e-9; %1360                                              %[m] center wavelength
        s.f0=const.c/s.wavelength;                                          %[1/s] center frequency
        s.w0=2*pi*s.f0;
        s.n_cycles=3;     %0.3657                                           % number of cycles inside pulse duration
        s.t_pulse=s.n_cycles/s.f0;                                          %[s] Pulse duration @ Intensity/e2
        s.t_fwhm= s.t_pulse*sqrt(log(2)/2);                                 %[s] Pulse duration @ Full-Width-Half-Maximum
        %Note: Scaling for const Peak Intensity for any Cyclenumber
        s.Q_In=0.4e-3*(s.n_cycles);                                      %[J] Pulse energy going in %1.9085e-3/3 for 14.79% ionization% 1.1776e-3/ncycle @1/2
        s.Q_Out=1e-3;                                                       %[J] Pulse energy at end of propagation        
        s.r_mode=100e-6; %30                                                   %[m] beam radius
        s.area_mode=pi*(s.r_mode)^2;                                        %[m^2] beam area
        s.Fluence=s.Q_In/(s.area_mode/2);                                   %[J/m^2] as peak Fluence is double for Gaussian shape!!!
        s.alpha=log(s.Q_In/s.Q_Out);                                        % Attenuation coefficient - Eloss
        end
    end
end