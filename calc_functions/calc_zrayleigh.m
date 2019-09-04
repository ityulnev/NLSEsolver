%This class calculates the rayleigh length for given pulse width
classdef calc_zrayleigh
    
    properties
    win, zr, Zsteps
    end
    
   methods 
       function s=calc_zrayleigh(beam,mesh,medium,pulse,waist)
        %initial beam waist at Intensity/e^2
        if waist>0
            s.win=waist;             
        else
            indexr=find(medium.Iconst.*abs(pulse.Ert(:,mesh.tmid)).^2<max(medium.Iconst.*abs(pulse.Ert(:,mesh.tmid)).^2)./exp(2),1);
            s.win=indexr.*mesh.dr+mesh.dr*2;
        end
        %Rayleigh length
            s.zr=pi*s.win^2/(beam.wavelength);
        s.Zsteps=round(s.zr/mesh.dz);         
       end
   end
end

    
      