
pp3=zeros(1,30);
pp2=0;
figure;
for hh=35900:5:35930

pp=polyfit(mesh.r(2:30)',angle(Erf(2:30,hh)),3);
pp2=polyval(pp,mesh.r(1:30));
pp3=[pp3;pp2];
hold on
plot(mesh.r(1:30),pp2)
pause(0.1)
end
figure;
plot(mesh.f(35900:5:35930),[pp3(2:end,1)])

hold on; plot(mesh.f(35900:5:35930),angle(Erf(1:4,35900:5:35930)))


