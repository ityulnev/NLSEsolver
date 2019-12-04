%Create structure that contains properties of a pulse in t,f and r domain
%Pulse made as sum of 2 electric fields in time domain
classdef sumpulse_init
    
    properties
    wavelength,f0,w0,t_delay,beam_area,Iconst,Ert,Energy,Erf,Energyf,pfmid,ptmid,fwhmF,fwhmT,Ipeak,t_Ie2,indt_Ie2
    end
    
    methods
        function s=sumpulse_init(mesh,p1,p2,Iconst,beam_area)                   
        %% pulse Field and Intensity in t,f
        s.Ert=p1.Ert+p2.Ert;
        s.Erf=myfft(s.Ert,mesh);
        %% Intensity
        s.Iconst=Iconst;
        Irt=s.Iconst.*abs(s.Ert).^2;
        Irf=s.Iconst.*abs(s.Erf).^2;         
        %find peak position in time and frequency
        s.ptmid=find(max(Irt(1,:))==Irt(1,:));
        s.pfmid=find(max(Irf(1,:))==Irf(1,:));
        %% Pulse Parameters
        s.t_delay=p1.t_delay+p2.t_delay;
%         p1.Energy+p2.Energy;
        s.f0=calc_centerofmass(mesh.f',Irf(1,:)','cartesian');
        s.w0=2*pi*s.f0;
        s.wavelength=const.c/s.f0;        
        %% Calculate Pulse energy  
        s.beam_area=beam_area;
        if size(s.Ert,1)>1
            s.Energy=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t,Irt,2),1);
            s.Energyf=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,Irf,2),1);
        else
            s.Energy=s.beam_area.*trapz(mesh.t,Irt); 
            s.Energyf=s.beam_area.*trapz(mesh.f,Irf); 
        end
        %% Peak Intensity
        [s.Ipeak,Ipeakpos]=max(max(Irt,[],1),[],2);
        %% Duration
        mybounds=find_bounds2(s.Ert(1,:));
        s.t_Ie2=mesh.dt*(mybounds(1,3)-mybounds(1,1)).*3;
        s.indt_Ie2=mybounds(1,3);
        s.fwhmF=calc_fwhm(mesh.f,Irf(1,:));
        s.fwhmT=calc_fwhm(mesh.t,Irt(1,:));
        end
    end
end