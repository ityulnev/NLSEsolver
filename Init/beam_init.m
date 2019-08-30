%Create structure with all needed parameters specific to the beam/pulse
classdef beam_init

    properties
    wavelength,f0,w0,Q_In,Q_Out,t_fwhm,r_mode,area_mode,Fluence,alpha,n_cycles
    end
    
    methods 
        function s=beam_init
        s.wavelength=2000e-9;%[m]
        s.f0=const.c/s.wavelength;
        s.w0=2*pi*s.f0;
        s.n_cycles=3;
        s.t_fwhm=s.n_cycles/s.f0;%[s] Pulse duration
        
        s.Q_In=2.1e-3*(s.n_cycles/3);%[J] Pulse energy going in
        %Note: Scaling for const Peak Intensity ~4.033e18W/m^2 for any Cyclenumber!
        s.Q_Out=1e-3;%[J] Pulse energy out        
        
%         factor=0.72;% some factor to vary the beams size - arbitrary!
%         trproduct=1.872e-18*sqrt(s.n_cycles/3);% Scaling for const peak Intensity leading to 1.24% ionization @ 3 cycles
        s.r_mode=100e-6;%trproduct/s.t_fwhm;%factor*(130e-6);%[m] beam radius
        s.area_mode=pi*(s.r_mode)^2;%[m^2] beam area
        s.Fluence=s.Q_In/(s.area_mode/2);%as peak Fluence is double for Gaussian shape!!!
        s.alpha=log(s.Q_In/s.Q_Out);% Attenuation coefficient - Eloss
        end
    end
end