%Structure holds finite difference matrix for (d/dx^2+d/dy^2) in cylindric
%coordinates which is split into vectors r(pos, neg, mid)
%Function: 1/r*(d/dr(r*d/dr(A)))=rneg*A(j-1)-(rneg+rpos)*A(j)+rpos*A(j+1)
function M_fd=M_FiniteDifference_init(mesh,boundary)
M_fd=struct;

M_fd.r_d=mesh.r(1:end-2);%r of regular length r0...r_Nr-2
M_fd.r_o=mesh.r(2:end-1);%r shifted up r1...r_Nr
M_fd.r_u=mesh.r(3:end);%r shifted up 2 times r2...r_Nr+1
M_fd.rneg=(M_fd.r_d+M_fd.r_o)./(2.*M_fd.r_o.*(mesh.dr^2));
M_fd.rpos=(M_fd.r_o+M_fd.r_u)./(2.*M_fd.r_o.*(mesh.dr^2));
% M_fd.rmid=-(M_fd.rneg+M_fd.rpos);
M_fd.boundary=boundary;
end