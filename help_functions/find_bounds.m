%Find left and right index at |E|/exp(1)
function [LRbound]=find_bounds(E)

LRbound=zeros(size(E,1),2);
for m=1:size(E,1)
Emax=max(abs(E(m,:).^2)./exp(2),[],2);
LRbound(m,1)=find(abs(E(m,:).^2)>Emax,1);
LRbound(m,2)=find(abs(E(m,:).^2)>Emax,1,'last');
end

end