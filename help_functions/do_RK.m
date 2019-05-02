%calculate next step in z via Runge Kutta
function [Erf]=do_RK(mesh,pulse,beam,medium,Erf,M_fd,m)
h=mesh.dz;
       k1= h.*calcfunctionRK(mesh,pulse,beam,medium,Erf,M_fd);   
       k2 = h.*calcfunctionRK(mesh,pulse,beam,medium,Erf+k1./2,M_fd);
       k3 = h.*calcfunctionRK(mesh,pulse,beam,medium,Erf+k2./2,M_fd);  
       k4 = h.*calcfunctionRK(mesh,pulse,beam,medium,Erf+k3,M_fd);       
Erf=Erf + (k1+2.*k2+2.*k3+k4)./6;
%%look at the 4 RK coefficients
% mid=pulse.pfmid; plot(mesh.r,[abs(k1(:,mid)),abs(k2(:,mid)),abs(k3(:,mid)),abs(k4(:,mid))])
% legend('k1','k2','k3','k4');
% title(['z=',num2str(m*mesh.dz)]);
% pause(0.1);
end