%Create structure with all needed parameters specific to the beam/pulse
classdef beam_init

    properties
    wavelength,f0,w0,Q_In,Q_Out,t_fwhm,r_mode,area_mode,Fluence,alpha
    end
    
    methods 
        function s=beam_init
        s.wavelength=800e-9;%[m]
        s.f0=const.c/s.wavelength;
        s.w0=2*pi*s.f0;
        s.Q_In=2.1e-3;%[J] Pulse energy in
        s.Q_Out=1e-3;%[J] Pulse energy out
        s.t_fwhm=35e-15;%[s] Pulse duration

        factor=0.5;% some factor to vary the beams size - arbitrary!
        s.r_mode=factor*(130e-6);%[m] beam radius
        s.area_mode=pi*(s.r_mode)^2;%[m^2] beam area
        s.Fluence=s.Q_In/s.area_mode;
        s.alpha=log(s.Q_In/s.Q_Out);% Attenuation coefficient - Eloss
        end
    end
end