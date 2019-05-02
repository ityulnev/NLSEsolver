% Implemented 4th order Low Storage Runge Kutta Algorithm
% X_j=X_j-1 + b_j * q_j
% with q_j=a_j * q_j-1 + h* f(x_j-1)
% Source: 
% (J. H. WILLIAMSON, 'Low-Storage Runge-Kutta Schemes',JOURNAL OF
% COMPUTATIONAL PHYSICS, 1980)
function [Erf]=do_LowStoreRK(mesh,pulse,beam,medium,Erf,M_fd,h)
%RK Coefficients
a_rk=[0, -5/9,-153/128];
b_rk=[1/3, 15/16, 8/15];
q=0;
%Calc Phase of reflected wave e(ikh)
    for mm=1:3
    q=h.*calc_mainfunctionRK(mesh,pulse,beam,medium,Erf,M_fd)+a_rk(mm).*q;   
    Erf=Erf+b_rk(mm).*q;
    [Erf]=set_boundaries(mesh,pulse,medium,Erf,M_fd,mm);
    end  
Erf(isnan(Erf))=0;
Erf(isinf(Erf))=0;
end




