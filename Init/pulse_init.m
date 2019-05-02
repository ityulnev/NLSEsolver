%Create structure that contains properties of the pulse in t,f and r domain
classdef pulse_init
    
    properties
    w0,t_fwhm,tau0,t0,Ert,carrier,A0,Energy,Erf,Energyf,order,Irf,Irt,pfmid,ptmid
    end
    
    methods
        function s=pulse_init(mesh,beam,medium,t0,order)
        s.order=order;                                                      %order: harmonic number
        s.t_fwhm=beam.t_fwhm/sqrt(s.order);                                 %Pulse Duration @ FWHM
        s.tau0=s.t_fwhm/(2*sqrt(log(2)));                                   %pulse duration @ tau0 (Gaussian variance)
        s.w0=beam.w0.*s.order;                                              %center frequency
        s.t0=t0;                                                            %time delay in s
        timedelay=-(t0*1i*2*pi.*(mesh.f));                                  %phase from time delay     
        s.carrier=1i*s.w0.*mesh.t;                                       %Oscillation of carrier wave with w0 center frequency
        %% calculate pulse
        ef=exp(-(2*pi.*(mesh.f)).^2.*s.tau0^2./2-timedelay);
        et=myifft(ef,mesh);
%         et=et.*exp(s.carrier);                                              %Excluded carrier wave for now...
        s.A0=sqrt(beam.Fluence/(sum(medium.Iconst.*abs(et).^2)*mesh.dt));
        n_gaussian=1;                                                       %Gaussian order
        er=exp(-(((mesh.r).^2./((beam.r_mode)^2))).^n_gaussian);
        %% pulse Field and Intensity in t,f
        s.Ert=s.A0.*transpose(er).*et;
        s.Erf=myfft(s.Ert,mesh);
        s.Erf=abs(s.Erf);
        s.Irt=medium.Iconst.*abs(s.Ert).^2;
        s.Irf=medium.Iconst.*abs(s.Erf).^2;
        
        %find peak position in time and frequency
        s.ptmid=find(max(s.Irt(1,:))==s.Irt(1,:));
        s.pfmid=find(max(s.Irf(1,:))==s.Irf(1,:));

        %% Calculate Pulse energy and Test       
        if size(s.Ert,1)>1
            s.Energy=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t,s.Irt,2),1);
            s.Energyf=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,s.Irf,2),1);
        else
            s.Energy=medium.area_hcf.*trapz(mesh.t,s.Irt); 
            s.Energyf=medium.area_hcf.*trapz(mesh.f,s.Irf); 
        end

        %test the energy conservation
        tolerance=1e-4;%set arbitrary tolerance
        test_errorMSG(abs(s.Energyf-s.Energy)/s.Energy >tolerance,'pulse_init: Energy of Et and Ef not conserved!')

        %test the time bandwidth product
        fwhmF=calc_fwhm(mesh.f,s.Irf(1,:));
        fwhmT=calc_fwhm(mesh.t,s.Irt(1,:));
        tolerance=2e-2;%set arbitrary tolerance
        test_errorMSG(abs(fwhmT*fwhmF-0.44)/0.44 >tolerance,'pulse_init: FWHM of Et and Ef not conserved!')
        end
    end
end