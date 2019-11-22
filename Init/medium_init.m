%Create structure that contains all properties of the medium in which the
%pulse propagates, e.g. gastype, refractive index, temperature, pressure
classdef medium_init
  
    properties
    temperature,pressure,n2,n0,n,Iconst,gas,k,k0,n0pressure,npressure,n_ext,npressure_ext,kGV,kGVD,kTOD,k1_w0,k2_w0,l,m,Eg,n_gas,P_crit,k_fit
    end
   
    methods
        function s=medium_init(mesh,beam,gas)
        s.gas=gas;%gastype
        s.temperature=270;%[K]
        s.pressure=1;%[bar]
        %% Refractive Index
        switch gas
            case 'Neon'
            [s.n0pressure,s.n0]=calc_refrIndex(beam.wavelength,gas,s.pressure,s.temperature);   %refractive index at center frequency
            [s.npressure,s.n]=calc_refrIndex(mesh.wvl,gas,s.pressure,s.temperature);            %refractive index for all frequencies
            %n2 Nonlinear Refractive Index pressure scaling
%             s.n2=s.pressure*0.85e-24;%0.625e-24;%[m^2/W]                                        %Nonlinear refractive index for SPM
            [n2pressure,n2]=calc_refrIndex(beam.wavelength,'Neon_n2',s.pressure,s.temperature);
            s.n2=(n2pressure-s.n0pressure)/2.224e18;
            case 'Argon'
            [s.n0pressure,s.n0] = calc_refrIndex(beam.wavelength,gas,s.pressure,s.temperature);   %refractive index at center frequency
            [s.npressure,s.n]   = calc_refrIndex(mesh.wvl,gas,s.pressure,s.temperature);    
            s.n2 = 9.5e-24; %[m^2/W]
        end
        %% Gas Parameters
        switch gas
            case 'Neon'
                    s.l=1;                                                           %Quantum Numbers l and m
                    s.m=0;
                    s.Eg=21.565*const.e;                                             %[J] from http://www.periodensystem.info/elemente/neon/
                    s.n_gas=2.686e25*s.pressure;                                                %1/m^3 
            case 'Argon'
                    s.l=1;
                    s.m=0;
                    s.Eg=15.76*const.e;                                              %[J] from http://www.periodensystem.info
                    s.n_gas=2.7e25*s.pressure;                                                  %1/m^3   
            case 'Xenon'
                    s.l=1;
                    s.m=0;
                    s.Eg=12.13*const.e;                                              %[J] from http://www.periodensystem.info
                    s.n_gas=2.4e25*s.pressure;                                                  %1/m^3  
        end
        %% Linear Refractive index and wave number k
        s.n_ext=[zeros(1,mesh.fbound-1),s.n];
        s.npressure_ext=[zeros(1,mesh.fbound-1),s.npressure];
        s.Iconst=0.5.*s.n0*const.c*const.eps0;% Constant factor for calculating signal intensity I=Iconst*abs(E)^2 ##0.5*
        % wave number k
        s.k0=s.n0pressure*beam.f0*2*pi/const.c;
        s.k=(s.npressure_ext).*(2.*pi.*mesh.f)./const.c;
        % expanded k=k0+k'+k''+k'''  
        s.kGV=gradient(s.k,mesh.df*2*pi);
        s.kGV(1:mesh.indexfmid)=0;
        s.kGVD=gradient(s.kGV,mesh.df*2*pi);
        s.kGVD(1:mesh.indexfmid+1)=0;
        %         s.kTOD=[diff(s.kGVD)./(mesh.df*2*pi),0];
        s.k1_w0=s.kGV(find(mesh.f>beam.f0,1)-1);
        s.k2_w0=s.kGVD(find(mesh.f>beam.f0,1)-1);
        %% Critical Power for Selffocusing>Divergence:
        alpha=1.8962;
        s.P_crit=alpha.*beam.wavelength^2/(4*pi*s.n0*s.n2);   
        %% k with reduced noise from linear fit
%         Vf=mesh.indexfmid:mesh.flength;
%         [aa]=polyfit(mesh.f(Vf),s.k(Vf),1);
%         s.k_fit=[zeros(1,mesh.indexfmid-1),mesh.f(Vf).*aa(1,1)];
%         
        
        end
    end
end