%Calculates FWHM of data y1
function [fwhmx]=calc_fwhm(x1,y1)
% Find the half max value.
halfMax = (min(y1) + max(y1)) / 2;
% Find where the data first drops below half the max.
index1 = find(y1 >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(y1 >= halfMax, 1, 'last');

fwhmx = x1(index2) - x1(index1);
end