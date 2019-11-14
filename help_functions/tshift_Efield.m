%Shifts electric field in time by dt
%E(t+dt)=A(t+dt) exp(i w0 [t+dt])
function [Et,Et2]=tshift_Efield(mesh,Et,dt)

Ef=myfft(Et,mesh);
Ef=Ef.*exp((1i*2*pi.*(mesh.f).*dt));
Et=myifft(Ef,mesh);

end