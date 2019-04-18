% Implemented 4th order Low Storage Runge Kutta Algorithm
% X_j=X_j-1 + b_j * q_j
% with q_j=a_j * q_j-1 + h* f(x_j-1)
% Source: 
% (J. H. WILLIAMSON, 'Low-Storage Runge-Kutta Schemes',JOURNAL OF
% COMPUTATIONAL PHYSICS, 1980)
function [Erf]=do_LowStoreRK(mesh,pulse,beam,fiber,Erf,M_fd,m)
%RK Coefficients
a_rk=[0, -5/9,-153/128];
b_rk=[1/3, 15/16, 8/15];
q=0;

Eprev=[Erf(1:3,:);Erf(end-2:end,:)];
%Calc Phase of reflected wave e(ikh)

    for mm=1:3
    q=mesh.dz.*calcfunctionRK(mesh,pulse,beam,fiber,Erf,M_fd)+a_rk(mm).*q;   
    Erf=Erf+b_rk(mm).*q;


    [Erf]=set_boundaries(mesh,fiber,Erf,M_fd,Eprev,mm);
%     fmid=find(abs(Erf(1,:))==max(abs(Erf(1,:))));
% kkr=log(Erf(2:end,fmid)./Erf(1:end-1,fmid))./(1i);
% % kkl=log(Erf(1:end-1,fmid)./Erf(2:end,fmid))./(1i);
% plot(mesh.r,[[0;real(kkr)./mesh.dr]]); legend([num2str(kkr(end)./mesh.dr)])
% xlim([0 140e-6])
% pause(0.1);

%     Eprev(4:6,:)=Erf(end-2:end,:);
    end


end




