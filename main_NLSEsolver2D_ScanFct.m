main_NLSEsolver2D_Fct_Scan
addpath('Init','help_functions','calc_functions','data')
% clear all
% close all
warning('none')                                                            %Reset previous warning messages
%% Read parameters from txt
% [param_var,param_val]=readfrom_txt('myparams.txt');
%[~,t_delay]=readfrom_txt('mydelay.txt');

%% initialize input Parameters
wavelength=2000e-9;
waist=100e-6;
Lz=1e-3;                                                                   %param_val(3);
pulse_cep=0*pi;                                                            %param_val(2)
pressure=5;                                                                %param_val(1);

mesh=mesh_init(wavelength,waist,Lz,2);  	                               % create meshgrid         
% pulse=general_pulse_init(mesh,wavelength,duration,waist,3.7e18,0,0);     % parameters of pulse      
% medium=medium_init(mesh,pulse,'Neon');                                   % parameters of medium   

%% Sum-Pulse
pulse1=general_pulse_init(mesh,1600e-9,7e-15,100e-6,2.1e18,0,t_delay,pulse_cep);
pulse2=general_pulse_init(mesh,800e-9,7e-15,100e-6,2.1e18,0,0,pulse_cep);
sumpulse=sumpulse_init(mesh,pulse1,pulse2,pulse1.Iconst,0); 
medium=medium_init_press(mesh,sumpulse,'Neon',pressure); 
% [n_e]=calc_2DeDensityADK(sumpulse.Ert(1,:),mesh,medium,sumpulse);
% IonizLvl=max(n_e(1,:))/medium.n_gas

%% 2D Propagate with Runge Kutta
comment='Synthesized pulse from 2 Gaussians';
boundcon=["open","const"]; %left and right boundary condition; 'open' or 'const'
[Er,Etrz,zprop,IonizLvl,Zsteps,mm,index_tL,index_tR]=do_WaveEqSolver(mesh,medium,sumpulse,sumpulse.Ert,boundcon,comment);

%% Save
if t_delay==0e-15
save([date,'sumpulse_scan',num2str(t_delay.*1e16),'dt_',num2str(param_val(2)),'piCEP_',num2str(medium.pressure),'bar.mat'],'mesh','medium','sumpulse','waist','Er','zprop','Zsteps','mm','-v7.3'); 
save([date,'sumpulse_scan',num2str(t_delay.*1e16),'dt_',num2str(param_val(2)),'piCEP_',num2str(medium.pressure),'bar_Etrz.mat'],'Etrz','index_tL','index_tR','-v7.3');
else
save([date,'sumpulse_scan',num2str(t_delay.*1e16),'dt_',num2str(param_val(2)),'piCEP_',num2str(medium.pressure),'bar.mat'],'sumpulse','waist','Er','zprop','Zsteps','mm','-v7.3');
save([date,'sumpulse_scan',num2str(t_delay.*1e16),'dt_',num2str(param_val(2)),'piCEP_',num2str(medium.pressure),'bar_Etrz.mat'],'Etrz','index_tL','index_tR','-v7.3');
end        