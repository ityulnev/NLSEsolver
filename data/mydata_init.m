%Structure contains data from a given matrix which is normalized and
%interpolated to the given mesh grid
classdef mydata_init
    
    properties
        E_IR,E_NIR,E_synth
    end
    
    methods
        function s=mydata_init(mesh,dataname)
           mydata=load(dataname);
    
           Enames = {'E_IR','E_NIR','E_synth'};
           fnames = fieldnames(mydata);
           for m=1:numel(fnames)-1
               if( isnumeric(mydata.(fnames{m})) )
                    s.(Enames{m}) = interp1(mydata.t_fs.*1e-15,mydata.(fnames{m}),mesh.t)./max(mydata.(fnames{m}));
               end
           end
                          
        end
    end
end
