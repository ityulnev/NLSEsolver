addpath('Init','help_functions','calc_functions','data')
% clear all
% close all
warning('none')% reset potential previous warnings
%% Read parameters from txt
[param_var,param_val]=readfrom_txt('myparams.txt');
[~,t_delay]=readfrom_txt('mydelay.txt');
%% initialize all input Parameters
wavelength=2000e-9;
waist=100e-6;
Lz=param_val(3);
pressure=param_val(1);
mesh=mesh_init(wavelength,waist,Lz,2);                                   % create meshgrid          #Set(r,t,f) #In(wavelength,beam_waist,Length in z,N Dimensions 1 or 2)   
% pulse=general_pulse_init(mesh,wavelength,duration,waist,3.7e18,0,0);     % parameters of pulse          #Set(Et(t,r,z=0),Ef,It,If)
% medium=medium_init(mesh,pulse,'Neon');                                   % parameters of medium     # Set(pressure,Temperature,Gas) -> Refractive index n(f),k,k1,k2...

%% 2 Gaussians electric field
% delay=(0).*1e-15;
% figure; hold on;
% for m=1:length(delay)
% pulse1=general_pulse_init(mesh,1600e-9,7e-15,100e-6,3.7e18,0,0);             %IR
% pulse2=general_pulse_init(mesh,800e-9,7e-15,100e-6,2.5e18,0,delay(m));     %NIR
% sumpulse=sumpulse_init(mesh,pulse1,pulse2,pulse1.Iconst,0); 
% [n_e]=calc_2DeDensityADK(sumpulse.Ert(1,:),mesh,medium,beam,sumpulse);
% IonizLvl=max(n_e(1,:))/medium.n_gas;
% plot(mesh.t,n_e(1,:)./medium.n_gas);
% end
% legend('0fs','2fs','4fs','6fs','8fs','10fs')

%% Sum Pulse
t_delay=1e-15;
my_CEP=param_val(2)*pi;
pulse1=general_pulse_init(mesh,1600e-9,7e-15,100e-6,2.1e18,0,t_delay,my_CEP);
pulse2=general_pulse_init(mesh,800e-9,7e-15,100e-6,2.1e18,0,0,my_CEP);
% figure; plot(mesh.t,[real(pulse2.Ert(1,:));abs(pulse2.Ert(1,:))]); xlim([-2e-14,2e-14])
sumpulse=sumpulse_init(mesh,pulse1,pulse2,pulse1.Iconst,0); 
medium=medium_init_press(mesh,sumpulse,'Neon',pressure); 


[n_e]=calc_2DeDensityADK(sumpulse.Ert(1,:),mesh,medium,sumpulse);
% plot(mesh.t,[n_e(1,:)./medium.n_gas]);
%figure; plot(mesh.t,real([pulse1.Ert(1,:);pulse2.Ert(1,:)]).^2)
IonizLvl=max(n_e(1,:))/medium.n_gas

%% 2D Propagate, Finite Difference + Split Step
comment='electron drift';
boundcon=["open","const"]; %left and right boundary condition
[Er,Etrz,zprop,IonizLvl,Zsteps,mm,index_tL,index_tR]=do_WaveEqSolver(mesh,medium,sumpulse,sumpulse.Ert,boundcon,comment);
%% Save
save([date,'sumpulse','_inmatlab.mat'],'-v7.3');
                    