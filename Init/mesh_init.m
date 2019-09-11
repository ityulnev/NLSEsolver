%awesome constructor for my mesh
classdef mesh_init
    
    properties
    Lz,dz,z,zlength,fmax,df,indexfmid,f,flength,fmin,freal,f0bound,t,dt,R,dr,r,rlength,wvl,dwvl,rmin,fbound,rmid,tmid,test,mode,tanhfilterL,tanhfilterR,tanhfilterLR
    end
    
   methods 
       function s=mesh_init(beam,Lz,dim)
        %z propagation
        s.Lz=Lz;%[m]
        s.dz=0.1e-6;
        s.z=0:s.dz:s.Lz;
        s.zlength=length(s.z);
        %frequency domain
        s.fmax=beam.f0*150;%[1/s]
        s.df=4e11;
        s.fmin=1e11;
        s.f=-s.fmax:s.df:s.fmax;
        s.flength=length(s.f);
        s.indexfmid=round(s.flength/2);
        s.freal=s.f+beam.f0;
        s.fbound=find(s.f>0,1,'first');
        s.f0bound=find(s.f>beam.f0,1,'first')-1;
        %wavelength
        s.wvl=const.c./s.f(s.fbound:end);
        s.dwvl=abs(s.wvl(3)-s.wvl(2));
        %time domain
        s.t=linspace(-1/(2*s.df),1/(2*s.df),s.flength);%[s]
        s.dt=abs(s.t(2)-s.t(1));
        s.tmid=round(length(s.t)/2);
        %radial / transverse 
        switch dim
            case 1
                s.r=0;
            case 2%3D with cylinder symmetry!
                s.R=150e-6;%800e-6;%[m]
                s.dr=2e-6;
                s.rmin=s.dr*2;
                s.r=s.rmin:s.dr:s.R;%start at r0=3*dr to avoid singularity at r0=0! 
        end
        if s.R<beam.r_mode
                %warning('pulse_init: FWHM of Et and Ef not conserved!') 
                dd=errordlg('mesh_init: !','Warning');
                uiwait(dd) 
        end
        s.rlength=length(s.r); 
        s.rmid=round(s.rlength/2);
       %% tanh(x) filter
       s.tanhfilterL=calc_tanhfilter((2*pi.*s.f-beam.w0).*s.t,s.indexfmid);
       s.tanhfilterR=fliplr(ifftshift(s.tanhfilterL));
       s.tanhfilterLR=s.tanhfilterL.*s.tanhfilterR;
       end
   end
end