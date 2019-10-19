%Calculates center of mass of data y1 in coordinate system of choice
function [com]=calc_centerofmass(x1,y1,csystem)

if isequal(size(x1), size(y1))
    %nothing
elseif (isvector(x1) && isvector(y1) && numel(x1) == numel(y1))
    x1=x1';
else
    warning('Error in calc_centerofmass(): x and y arrays dont match');
end
    
switch csystem
    case 'cylinder'    
        com = sqrt(trapz(x1,x1.^3.*y1,1)./trapz(x1,x1.*y1,1));
    case 'cartesian'    
        com = trapz(x1,x1.*y1,1)./trapz(x1,y1,1);      
end
end