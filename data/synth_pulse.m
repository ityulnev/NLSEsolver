%Create structure that contains properties of the pulse in t,f and r domain
classdef synth_pulse
    
    properties
       Ert,Erf,E_IR,E_NIR,E_NIR_init,f0,gaussfilter,t_gfilter,Iconst,I_peak,t_delay,r_mode1,r_mode2,Energy,Energyf,ptmid,pfmid
    end
    
    methods
        function s=synth_pulse(mesh,E1,E2,r_mode1,r_mode2,I_peak,t_delay,t_gfilter)

        %% pulse Field and Intensity in t,f
        s.t_gfilter=t_gfilter;
        s.gaussfilter=calc_supergaussian(mesh.t,s.t_gfilter,10,0); 
        Et1=get_compEField(mesh,E1).*s.gaussfilter;
        Et2=get_compEField(mesh,E2).*s.gaussfilter;
        E_IR=do_filter(Et1,'tanhfilterLR','inF',mesh);
        E_NIR=do_filter(Et2,'tanhfilterLR','inF',mesh);
        %Field ratios from I_1=3 * I_2
        s.Iconst=const.eps0*const.c*1/2;
        s.I_peak=I_peak;
        factor_ItoE=s.I_peak/s.Iconst;
        s.E_IR=sqrt(factor_ItoE)*(sqrt(3)/(1+sqrt(3))).*E_IR;
        s.E_NIR_init=sqrt(factor_ItoE)*(1/(1+sqrt(3))).*E_NIR;
        %% Time shift
        s.t_delay=t_delay;
        s.E_NIR=tshift_Efield(mesh,s.E_NIR_init,s.t_delay);
        %% Radial component                        
        s.r_mode1 = r_mode1;
        s.r_mode2 = r_mode2;
        er1 = exp(-(((mesh.r).^2./((s.r_mode1)^2))));
        er2 = exp(-(((mesh.r).^2./((s.r_mode2)^2))));
        %IR
        s.Ert = transpose(er1).*s.E_IR+transpose(er2).*s.E_NIR;
        s.Erf = myfft(s.Ert,mesh);

        %% Intensity
        Irt=s.Iconst.*abs(s.Ert).^2;
        Irf=s.Iconst.*abs(s.Erf).^2;         
        %find peak position in time and frequency
        s.ptmid=find(max(Irt(1,:))==Irt(1,:));
        s.pfmid=find(max(Irf(1,:))==Irf(1,:));
        %
        s.f0=(trapz(mesh.f,Irf(1,:).*mesh.f)./trapz(mesh.f,Irf(1,:)));

        %% Calculate Pulse energy and Test       
        if size(s.Ert,1)>1
            s.Energy=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t,Irt,2),1);
            s.Energyf=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,Irf,2),1);
        else
            s.Energy=pi*s.r_mode1^2.*trapz(mesh.t,Irt); 
            s.Energyf=pi*s.r_mode2^2.*trapz(mesh.f,Irf); 
        end
        %test the energy conservation
        tolerance=1e-4;%set arbitrary tolerance
        energy_error=abs(s.Energyf-s.Energy)/s.Energy;
%         test_errorMSG(fwhm_error>tolerance,'pulse_init: Energy of Et and Ef not conserved!')
        if energy_error>tolerance
            warning(['pulse_init: Energy conservation between Ef and Et deviate by ',num2str(energy_error.*100),'%'])    
        end
        end
    end
end