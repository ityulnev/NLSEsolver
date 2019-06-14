%Calculates a Gaussian function in xdata 
%With variance/2 @E^2/exp(2)=halfwidth
%Super-Gaussian order N
function E=calc_supergaussian(xdata,halfwidth,order)
E=exp(-(((xdata).^2./((halfwidth).^2))).^order);
end