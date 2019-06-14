%Test for given condition
%for condition==true show error message and stop code
function []= test_errorMSG(condition,message)  
    if condition
        dd=errordlg(message,'Warning');
        uiwait(dd)
        warning('myErrorID')   
    end
end