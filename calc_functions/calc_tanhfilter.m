%Calculates one sided tanh(x) filter 
%[0:x0] -> tanh(x)=0
%[x0+1:tanhmax-1] -> tanh(x)
%[tanhmax:end] -> tanh(x)=1
function [tanhfilter]= calc_tanhfilter(x1,x0)
tanhfilter=tanh(x1);
tanhmax=find(max(tanhfilter)==tanhfilter,1);
tanhfilter(tanhmax:end)=1;
tanhfilter(1:x0)=0;
end