%This class calculates pulse width Intensity/e^2 for given pulse width
%after z distance of propagation
function [wprop]=calc_wprop(waist,zr,zlength)
    wprop=waist*sqrt(1+(zlength/zr)^2);
end

