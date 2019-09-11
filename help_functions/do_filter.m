%Applies filter function of choice onto matrix/field
%domain=='inF' applied in frequency domain
%else regular filter
function [Er]=do_filter(Er,filtertype,domain,mesh)

switch filtertype
    case 'tanhfilterLR'
        myfilter=mesh.tanhfilterLR;
    case 'tanhfilterL'
        myfilter=mesh.tanhfilterL;
    case 'tanhfilterR'
        myfilter=mesh.tanhfilterR;       
end
if strcmp('inF',domain)
    Er=myifft(myfilter.*myfft(Er,mesh),mesh);
else
    Er=myfilter.*Er;
end

end