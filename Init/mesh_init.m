%Construct a mesh-grid in dimensions (t<->f,r,z) 
%t-time
%f-frequency
%r-radial (Cyllinder symmetrie)
%z-propagation direction
classdef mesh_init
    
    properties
    wavelength,f0,Lz,dz,z,zlength,fmax_fac,fmax,df,indexfmid,f,flength,fmin,freal,f0bound,t,dt,R,dr,r,rlength,wvl,dwvl,rmin,fbound,rmid,tmid,Tfilter_fac,Tfilter_L,Tfilter_R,Tfilter_LR,Gfilter_Thalfwidth,Gfilter_T,Gfilter_Rhalfwidth,Gfilter_R,Gfilter_fac
    end
    
   methods 
       function s=mesh_init(wavelength,r_mode,Lz,dim)
        %z propagation
        s.Lz=Lz;%[m]
        s.dz=0.5e-6;
        s.z=0:s.dz:s.Lz;
        s.zlength=length(s.z);
        %frequency domain
        s.wavelength=wavelength;
        s.f0=const.c/s.wavelength; %f mesh based on center frequency f0!
        s.fmax_fac=100;
        s.fmax=s.f0*s.fmax_fac;%[1/s]
        s.df=4e11;
        s.fmin=1e11;
        s.f=-s.fmax:s.df:s.fmax;
        s.flength=length(s.f);
        s.indexfmid=round(s.flength/2);
        s.freal=s.f+s.f0;
        s.fbound=find(s.f>0,1,'first');
        s.f0bound=find(s.f>s.f0,1,'first')-1;
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
                s.R=180e-6;%800e-6;%[m]
                s.dr=1.5e-6;
                s.rmin=s.dr*2;
                s.r=s.rmin:s.dr:s.R;%start at r0=3*dr to avoid singularity at r0=0! 
        end
        if s.R<r_mode
                dd=errordlg('mesh_init: Beam radius bigger than r-mesh!','Warning');
                uiwait(dd) 
        end
        s.rlength=length(s.r); 
        s.rmid=round(s.rlength/2);
       %% tanh(x) & Supergaussian(n=10) filters
       s.Tfilter_fac=4;
       s.Tfilter_L=calc_tanhfilter((s.f0.*s.t.*s.Tfilter_fac),s.indexfmid);
       s.Tfilter_R=fliplr(ifftshift(s.Tfilter_L));
       s.Tfilter_LR=s.Tfilter_L.*s.Tfilter_R;
       s.Gfilter_Thalfwidth=800e-15;
       s.Gfilter_T=calc_supergaussian(s.t,s.Gfilter_Thalfwidth,10,0);
       s.Gfilter_fac=1;
       s.Gfilter_Rhalfwidth=100e-6*2;
       s.Gfilter_R=calc_supergaussian(s.r,s.Gfilter_Rhalfwidth,10,0);
       end
   end
end