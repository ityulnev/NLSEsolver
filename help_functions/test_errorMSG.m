%Test for given condition
%for condition==true show error message and stop code
function []= test_errorMSG(condition,message)  
    if condition
        switch globproperties.mode 
            case 'debug'
                dd=errordlg(message,'Warning');
                uiwait(dd)
            case 'cluster'
                %nothing     
        end
        warning('myErrorID')   
    end
end