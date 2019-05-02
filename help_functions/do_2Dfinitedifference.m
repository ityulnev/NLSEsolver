% Calculates transverse derivative with Finite Difference Method in cylinder
% -coordinates
% [(d/dx)^2+ (d/dy)^2]A=1/r * d/dr[ r * d/dr(A)]
function [Erf]=do_2Dfinitedifference(mesh,medium,Erf,M_fd)         
% E_d=Erf(1:end-2,:); %A_0...A_Nr-1
% E_o=Erf(2:end-1,:); %A_1...A_Nr
% E_u=Erf(3:end,:); %A2...A_Nr+1
newErf=(M_fd.rneg'.*Erf(1:end-2,:)-(M_fd.rneg+M_fd.rpos)'.*Erf(2:end-1,:)+M_fd.rpos'.*Erf(3:end,:));
boundL=newErf(1,:);
boundR=newErf(end,:);
Erf=[boundL;newErf;boundR];
end