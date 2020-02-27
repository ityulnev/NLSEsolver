%My main function for two dimensional pulse propagation in gas 
%Used for parameter scan f(pulse-delay,carrier-envelope-phase)
%M-Shaped pulse spectrum built from two Gaussians at 800nm and 1600nm
function main_NLSEsolver2D_Scan(t_delay,cep)
addpath('Init','help_functions','calc_functions','data')
% clear all
% close all
warning('none')
%% Read parameters from txt
[param_var,param_val]=readfrom_txt('myparams.txt');
%[~,t_delay]=readfrom_txt('mydelay.txt');

%% initialize all input Parameters
wavelength=2000e-9;
waist=100e-6;
Lz=param_val(3);
pressure=param_val(1);
mesh=mesh_init(wavelength,waist,Lz,2);  	                               % create meshgrid          #Set(r,t,f) #In(wavelength,beam_waist,Length in z,N Dimensions 1 or 2)   
% pulse=general_pulse_init(mesh,wavelength,duration,waist,3.7e18,0,0);         % parameters of pulse          #Set(Et(t,r,z=0),Ef,It,If)
% medium=medium_init(mesh,pulse,'Neon');                                       % parameters of medium     # Set(pressure,Temperature,Gas) -> Refractive index n(f),k,k1,k2...
pulse_cep=cep*pi;

%% Sum Pulse
pulse1=general_pulse_init(mesh,1600e-9,7e-15,100e-6,2.1e18,0,t_delay,pulse_cep);
pulse2=general_pulse_init(mesh,800e-9,7e-15,100e-6,2.1e18,0,0,pulse_cep);
sumpulse=sumpulse_init(mesh,pulse1,pulse2,pulse1.Iconst,0); 
medium=medium_init_press(mesh,sumpulse,'Neon',pressure); 
% [n_e]=calc_2DeDensityADK(sumpulse.Ert(1,:),mesh,medium,sumpulse);
% IonizLvl=max(n_e(1,:))/medium.n_gas

%% 2D Propagate with Runge Kutta
comment='electron drift';
boundcon=["open","const"]; %left and right boundary condition
[Er,Etrz,zprop,IonizLvl,Zsteps,mm,index_tL,index_tR]=do_WaveEqSolver(mesh,medium,sumpulse,sumpulse.Ert,boundcon,comment);

%% Save
if t_delay==0e-15
save([date,'sumpulse_scan',num2str(t_delay.*1e16),'dt_',num2str(cep*10),'piCEP_',num2str(medium.pressure),'bar.mat'],'mesh','medium','sumpulse','waist','Er','zprop','Zsteps','mm','-v7.3'); 
save([date,'sumpulse_scan',num2str(t_delay.*1e16),'dt_',num2str(cep*10),'piCEP_',num2str(medium.pressure),'bar_Etrz.mat'],'Etrz','index_tL','index_tR','-v7.3');
else
save([date,'sumpulse_scan',num2str(t_delay.*1e16),'dt_',num2str(cep*10),'piCEP_',num2str(medium.pressure),'bar.mat'],'sumpulse','waist','Er','zprop','Zsteps','mm','-v7.3');
save([date,'sumpulse_scan',num2str(t_delay.*1e16),'dt_',num2str(cep*10),'piCEP_',num2str(medium.pressure),'bar_Etrz.mat'],'Etrz','index_tL','index_tR','-v7.3');
end        
            
 %% Save
%if t_delay==0e-15
%save([date,'sumpulse_scan',num2str(t_delay.*1e16),'dt_',num2str(param_val(2)),'piCEP_',num2str(medium.pressure),'bar.mat'],'mesh','medium','sumpulse','waist','Er','zprop','Zsteps','mm','Etrz','index_tL','index_tR','-v7.3'); 
%else
%save([date,'sumpulse_scan',num2str(t_delay.*1e16),'dt_',num2str(param_val(2)),'piCEP_',num2str(medium.pressure),'bar.mat'],'sumpulse','waist','Er','zprop','Zsteps','mm','Etrz','index_tL','index_tR','-v7.3');
%end   

end       