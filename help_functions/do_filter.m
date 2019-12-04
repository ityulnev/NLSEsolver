%Applies filter function of choice onto matrix/field
%domain=='inF' applied in frequency domain
%else regular filter
function [Er]=do_filter(Er,filtertype,domain,mesh)

switch filtertype
    case 'tanhfilterLR'
        myfilter=mesh.Tfilter_LR;
    case 'tanhfilterL'
        myfilter=mesh.Tfilter_L;
    case 'tanhfilterR'
        myfilter=mesh.Tfilter_R;  
    case 'tanhOuterLR'
        tempfilter=mesh.Tfilter_R;
        tempfilter(1:mesh.indexfmid)=1;
        myfilter=tempfilter.*fliplr(tempfilter);
end
if strcmp('inF',domain)
    Er=myifft(myfilter.*myfft(Er,mesh),mesh);
else
    Er=myfilter.*Er;
end

end