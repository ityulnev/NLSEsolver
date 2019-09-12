%Calculates Envelope from real electric field
%as Envelope=abs(Ecomplex)
function [E]=get_envelope(E,mesh)
E=abs(myifft(get_compfield(myfft(E,mesh),mesh),mesh));
end

function Ef=get_compfield(Ef,mesh)
Ef=Ef.*2;
Ef(1:mesh.indexfmid)=0;
end