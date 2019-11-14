%Take real electric field E(t) and turn it complex
%Real Field: Ereal = A exp(i w0 t) + c.c.
%Complex Field: Ecomp = A exp(i w0 t)
function [Et]=get_compEField(mesh,Et)
Ef = myfft(Et,mesh);
Ef(1:mesh.indexfmid)=0;
Ef = Ef.*2;
Et = myifft(Ef,mesh);
end
