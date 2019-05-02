%Create structure that contains all properties of the medium in which the
%pulse propagates, e.g. gastype, refractive index, temperature, pressure
classdef medium_init
  
    properties
    temperature,pressure,n2,n0,n,Iconst,gas,k,k0,n0pressure,npressure,n_ext,npressure_ext,kGV,kGVD,kTOD
    end
   
    methods
        function s=medium_init(mesh,beam,gas)
        s.gas=gas;%gastype
        s.temperature=300;%[K]
        s.pressure=2;%[bar]
        switch gas
            case 'Neon'
            s.n2=s.pressure*0.85e-24;%0.625e-24;%[m^2/W]                                        %Nonlinear refractive index for SPM
            [s.n0pressure,s.n0]=calc_refrIndex(beam.wavelength,gas,s.pressure,s.temperature);   %refractive index at center frequency
            [s.npressure,s.n]=calc_refrIndex(mesh.wvl,gas,s.pressure,s.temperature);            %refractive index for all frequencies
        end        
        %% nefractive index at negative frequencies = 0
        dfbound=find(mesh.f==beam.f0,1)-mesh.fbound;%shift f0 into the middle of array
        s.n_ext=[zeros(1,(mesh.fbound-1)-dfbound),s.n,zeros(1,dfbound)];
        s.npressure_ext=[zeros(1,(mesh.fbound-1)-dfbound),s.npressure,zeros(1,dfbound)];%[zeros(1,mesh.fbound-1),s.npressure];

        s.Iconst=0.5*s.n0*const.c*const.eps0;% Constant factor for calculating signal intensity I=Iconst*abs(E)^2
        % wave number k
        s.k=s.npressure_ext.*(2.*pi.*mesh.f)./const.c;
        % expanded k=k0+k'+k''+k'''  
        s.k0=s.n0pressure*beam.f0*2*pi/const.c;
        s.kGV=[0,diff(s.k)./(mesh.df*2*pi)]; 
        s.kGV((mesh.fbound)-dfbound)=0;
        s.kGV(mesh.flength-dfbound+1)=0;
        s.kGVD=[diff(s.kGV)./(mesh.df*2*pi),0];
        s.kGVD((mesh.fbound)-dfbound)=0;
        s.kGVD(mesh.flength-dfbound)=0;
%         s.kTOD=[diff(s.kGVD)./(mesh.df*2*pi),0];

        end
    end
end