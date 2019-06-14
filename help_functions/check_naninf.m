%Function look through NxN matrix and find number of Not a Number and
%Infinite values
function [Nnan,Ninf]=check_naninf(Matrix)
Nnan=sum(sum(isnan(Matrix),1),2);
Ninf=sum(sum(isinf(Matrix),1),2);
end