%Calculates a Gaussian function in xdata 
%With variance/2 @E^2/exp(2)=halfwidth
%Super-Gaussian order N
%offset is t0 from exp(-(t-t0/tau).^2.^N)
function E=calc_supergaussian(xdata,halfwidth,order,offset)
E=exp(-(((xdata-offset).^2./((halfwidth).^2))).^order);
end