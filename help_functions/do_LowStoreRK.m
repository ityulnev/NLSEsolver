% Implemented 4th order Low Storage Runge Kutta Algorithm
% X_j=X_j-1 + b_j * q_j
% with q_j=a_j * q_j-1 + h* f(x_j-1)
% Source: 
% (J. H. WILLIAMSON, 'Low-Storage Runge-Kutta Schemes',JOURNAL OF
% COMPUTATIONAL PHYSICS, 1980)
function [Er]=do_LowStoreRK(mesh,pulse,medium,Er,M_fd,h)
%RK Coefficients
a_rk=[0, -5/9,-153/128];
b_rk=[1/3, 15/16, 8/15];
q=0;
%Calc Phase of reflected wave e(ikh)
Eprev=[Er(1:2,:);Er(end-1:end,:)];
    for mm=1:3
    q=h.*calc_mainfunctionRK(mesh,pulse,medium,Er,M_fd)+a_rk(mm).*q;   
    Er=Er+b_rk(mm).*q;
    [Er]=set_boundaries(mesh,pulse,medium,Er,M_fd,mm,Eprev);
    end  
end




