%Create structure that contains properties of the pulse in t,f and r domain
classdef pulse_init
    
    properties
    w0,t_pulse,tau0,t0,Ert,carrier,A0,Energy,Erf,Energyf,order,pfmid,ptmid,fwhmF,fwhmT,Epeak,Ipeak,PpeakTheo,IpeakTheo
    end
    
    methods
        function s=pulse_init(mesh,beam,medium,t_delay,order)
        s.order=order;                                                     %order: harmonic number
        s.t_pulse=beam.t_pulse/sqrt(s.order);                              %Pulse duration @ I/e^2
        s.tau0=s.t_pulse/(2*sqrt(2));                                      %Pulse duration @ 1sigma (Gaussian variance)
        s.w0=beam.w0.*s.order;                                             %center frequency
        s.t0=t_delay;                                                      %time delay in s
        timedelay=-(t_delay*1i*2*pi.*(mesh.f));                            %phase from time delay     
        s.carrier=1i*s.w0.*(mesh.t-s.t0);                                  %Oscillation of carrier wave with w0 center frequency
        %% calculate pulse
        ef=exp(-(2*pi.*(mesh.f)).^2.*s.tau0^2./2-timedelay);
        et=myifft(ef,mesh);
        et=et.*exp(s.carrier);
        s.A0=sqrt(beam.Fluence/(sum(medium.Iconst.*abs(et).^2)*mesh.dt));
        n_gaussian=1;                                                      %Gaussian order
        er=exp(-(((mesh.r).^2./((beam.r_mode)^2))).^n_gaussian);
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
        Irt=medium.Iconst.*abs(s.Ert).^2;
        Irf=medium.Iconst.*abs(s.Erf).^2;         
        %find peak position in time and frequency
        s.ptmid=find(max(Irt(1,:))==Irt(1,:));
        s.pfmid=find(max(Irf(1,:))==Irf(1,:));
        %% Peak Intensity
        [s.Ipeak,Ipeakpos]=max(max(Irt,[],1),[],2);
        s.Epeak=max(abs(s.Ert(:,Ipeakpos)),[],1);
        s.PpeakTheo=0.94*(beam.Q_In/(s.t_pulse*sqrt(log(2)/2))/sqrt(order));
        s.IpeakTheo=s.PpeakTheo/(beam.area_mode/2);%Peak intensity for Gaussian beam from peak power
        %% Test Peak Intensity        
        tolerance=1e-2;
        Ipeak_error=abs(s.IpeakTheo-s.Ipeak)/s.IpeakTheo;
        if  Ipeak_error<tolerance
            %do nothing  
        else
            warning(['pulse_init: Peak Intensity of Irt deviates from theoretical Gaussian Ipeak by ',num2str(Ipeak_error.*100),'%'])    
        end        
        %% Calculate Pulse energy and Test       
        if size(s.Ert,1)>1
            s.Energy=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t,Irt,2),1);
            s.Energyf=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,Irf,2),1);
        else
            s.Energy=medium.area_hcf.*trapz(mesh.t,Irt); 
            s.Energyf=medium.area_hcf.*trapz(mesh.f,Irf); 
        end
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
        tolerance=2e-2;%set arbitrary tolerance
        fwhm_error=abs(s.fwhmT*s.fwhmF-0.44)/0.44;
%         test_errorMSG(abs(s.fwhmT*s.fwhmF-0.44)/0.44 >tolerance,'pulse_init: FWHM of Et and Ef not conserved!')  
        if fwhm_error>tolerance
%             warning(['pulse_init: FWHM product of Irt/Irf deviate from theoretical value by ',num2str(energy_error.*100),'%'])    
                    warning(['pulse_init: FWHM product of Irt/Irf deviate from theoretical value by ',num2str(fwhm_error.*100),'%'])    
        end
        end
    end
end