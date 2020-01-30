%Create structure that contains properties of the pulse in t,f and r domain
classdef general_pulse_init
    
    properties
    wavelength,f0,w0,t_pulse,tau0,t_delay,cep,Ipeak_in,r_mode,beam_area,Iconst,Ert,carrier,A0,Energy,Erf,Energyf,pfmid,ptmid,fwhmF,fwhmT,Epeak,Ipeak,PpeakTheo,IpeakTheo,t_Ie2,indt_Ie2
    end
    
    methods
        function s=general_pulse_init(mesh,wavelength,t_pulse,r_mode,Ipeak_in,Iconst,t_delay,cep)            
        s.t_pulse=t_pulse;                                                 %Pulse duration @ I/e^2
        s.tau0=s.t_pulse/(2*sqrt(2));                                      %Pulse duration @ 1sigma (Gaussian variance)
        s.wavelength=wavelength;
        s.f0=const.c/s.wavelength;
        s.w0=2*pi*s.f0;
        s.t_delay=t_delay;                                                 %time delay in s
        timedelay=-(s.t_delay*1i*2*pi.*(mesh.f));                          %phase from time delay     
        s.carrier=1i*s.w0.*(mesh.t-s.t_delay);                             %Oscillation of carrier wave with w0 center frequency
        s.cep=cep;
        if Iconst>0
            s.Iconst=Iconst;                                               %1/2*eps0*c*n0 factor for Calculating Intensity
       else
            s.Iconst=const.eps0*const.c*0.5;
       end
        %% calculate pulse
        ef=exp(-(2*pi.*(mesh.f)).^2.*s.tau0^2./2-timedelay);
        et=myifft(ef,mesh);
        et=et.*exp(1i.*cep).*exp(s.carrier)./max(abs(et));
        s.Ipeak_in=Ipeak_in;
        s.A0=sqrt(Ipeak_in/s.Iconst);  % integral over Envelope^2 = integral over time averaged Poynting vector!
                                                               % 0.5*int(abs(Ecomplex)^2,dt)dt=int(abs(real(Ecomplex))^2,dt)
        n_gaussian=1;                                                      %Gaussian order
        s.r_mode=r_mode;
        er=exp(-(((mesh.r).^2./((s.r_mode)^2))).^n_gaussian);
        %% pulse Field and Intensity in t,f
        s.Ert=s.A0.*transpose(er).*et;
        s.Erf=myfft(s.Ert,mesh);
        %% Handling of negative frequencies:
        %Gaussian mirror
        gmirror=-flip(s.Erf,2);
        s.Erf=s.Erf+gmirror; 
        s.Erf(:,1:mesh.indexfmid)=0;
        s.Ert=myifft(s.Erf,mesh);
        %% Intensity
        Irt=s.Iconst.*abs(s.Ert).^2;
        Irf=s.Iconst.*abs(s.Erf).^2;         
        %find peak position in time and frequency
        s.ptmid=find(max(Irt(1,:))==Irt(1,:));
        s.pfmid=find(max(Irf(1,:))==Irf(1,:));
        %% Energy
        s.beam_area=pi*s.r_mode^2;
        if size(s.Ert,1)>1
            s.Energy=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t,Irt,2),1);
            s.Energyf=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,Irf,2),1);
        else
            s.Energy=s.beam_area.*trapz(mesh.t,Irt); 
            s.Energyf=s.beam_area.*trapz(mesh.f,Irf); 
        end
        %% Peak Intensity
        [s.Ipeak,Ipeakpos]=max(max(Irt,[],1),[],2);
        s.Epeak=max(abs(s.Ert(:,Ipeakpos)),[],1);
        s.PpeakTheo=0.94*(s.Energy/(s.t_pulse*sqrt(log(2)/2)));%Peak Power of Gaussina pulse
        s.IpeakTheo=s.PpeakTheo/(s.beam_area/2);%Peak intensity for Gaussian beam from peak power
        %% Test Peak Intensity        
        tolerance=1e-2;
        Ipeak_error=abs(s.IpeakTheo-s.Ipeak)/s.IpeakTheo;
        if  Ipeak_error<tolerance
            %do nothing  
        else
            warning(['pulse_init: Peak Intensity of Irt deviates from theoretical Gaussian Ipeak by ',num2str(Ipeak_error.*100),'%'])    
        end        
        %% Test       
        %test the energy conservation
        tolerance=1e-4;%set arbitrary tolerance
        energy_error=abs(s.Energyf-s.Energy)/s.Energy;
%         test_errorMSG(fwhm_error>tolerance,'pulse_init: Energy of Et and Ef not conserved!')
        if energy_error>tolerance
            warning(['pulse_init: Energy conservation between Ef and Et deviate by ',num2str(energy_error.*100),'%'])    
        end
        %test the time bandwidth product
        s.fwhmF=calc_fwhm(mesh.f,Irf(1,:));
        s.fwhmT=calc_fwhm(mesh.t,Irt(1,:));
        mybounds=find_bounds2(s.Ert(1,:));
        s.t_Ie2=mesh.dt*(mybounds(1,3)-mybounds(1,1));
        s.indt_Ie2=mybounds(1,3);
        tolerance=2e-2;%set arbitrary tolerance
        fwhm_error=abs(s.fwhmT*s.fwhmF-0.44)/0.44;
%         test_errorMSG(abs(s.fwhmT*s.fwhmF-0.44)/0.44 >tolerance,'pulse_init: FWHM of Et and Ef not conserved!')  
        if fwhm_error>tolerance
%             warning(['pulse_init: FWHM product of Irt/Irf deviate from theoretical value by ',num2str(energy_error.*100),'%'])    
                    warning(['pulse_init: FWHM product of Irt/Irf deviate from theoretical value for a Gaussian Pulse by ',num2str(fwhm_error.*100),'%'])    
        end
        end
    end
end