%Find left and right index of general function for values left and right
function [LRbound]=find_mybounds(f,left,right)

LRbound=zeros(size(f,1),2);
for m=1:size(f,1)
LRbound(m,1)=find(f>left,1);
LRbound(m,2)=find(f<right,1,'last');
end

end