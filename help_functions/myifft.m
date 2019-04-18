%Energy conserving ifft (f->t)
function [Et]=myifft(Ef,mesh)
if size(Ef,1)>1
    Et=fftshift(ifft(ifftshift(Ef),[],2)).*(length(Ef)*mesh.df);
    if strcmp(mesh.test,'yes')
        intEt=trapz(mesh.r,trapz(mesh.t,abs(Et).^2,2),1);
        intEf=trapz(mesh.r,trapz(mesh.f,abs(Ef).^2,2),1);
    end
elseif size(Ef,1)==1
    Et=fftshift(ifft(ifftshift(Ef))).*(length(Ef)*mesh.df);  
    if strcmp(mesh.test,'yes')
        intEt=trapz(mesh.t,abs(Et).^2);
        intEf=trapz(mesh.f,abs(Ef).^2);
    end
else
    error('myifft: Dimension m of mxn matrix not real pos integer')
end
%Energy conservation check
if strcmp(mesh.test,'yes')
    tolerance=1e-4;
    test_errorMSG(abs(intEf-intEt)/intEf >tolerance,'myifft: Energy of Et and Ef not conserved!')
end
end



  