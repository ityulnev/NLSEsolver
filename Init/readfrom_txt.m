%Read parameters from txt file
function [param_var,param_val]=readfrom_txt(name)
paramtxt = fopen(name);
charInFile = textscan(paramtxt, '%s %n', 'delimiter', ',');
param_var=charInFile{1,1};
param_val=charInFile{1,2};
fclose(paramtxt);
end