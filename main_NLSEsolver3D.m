addpath('Init','help_functions','calc_functions','data')
% clear all
% close all
warning('none')
%% initialize all input Parameters
beam=beam_init;                        % parameters of the beam   #Set(wavelength, energy, duration)  
mesh=mesh_init(beam,1e-3,2);           % create meshgrid          #Set(r,t,f) #In(beam,propagate Distance Lz, Dimension t or t&r)   
medium=medium_init(mesh,beam,'Neon');  % parameters for medium    #Set(refractive index n, n2, beamsize)   
pulse=pulse_init(mesh,beam,medium,0,1,0);% parameters of pulse          #Set(Et(t,r,z=0),Ef,It,If)

% [n_e]=calc_2DeDensityADK(pulse.Ert(1,:),mesh,medium,beam,pulse);
% IonizLvl=max(n_e(1,:))/medium.n_gas

%% Measured electric fields
% dat=mydata_init(mesh,'Retrieved-May2019.mat');
% dt=0e-15;
% syn=synth_pulse(mesh,dat.E_IR,dat.E_NIR,30e-6,20e-6,2e18,dt,800e-15);
% figure; plot(mesh.t,[dat.E_synth.*3.876e10;real(syn.Ert(1,:))])
% figure; plot(mesh.f,abs([norm_fields(myfft(dat.E_synth,mesh),myfft(real(syn.Ert(1,:)),mesh),'indiv')]).^2)

%% 2D Propagate, Finite Difference + Split Step
comment='electron drift';
boundcon=["open","open"]; %left and right boundary condition
[Er,Etrz,zprop,IonizLvl,Zsteps,mm]=do_WaveEqSolver(mesh,beam,medium,pulse,pulse.Ert,boundcon,comment);
%% Self-focusing length
% L_spm=0.367*medium.k0*beam.r_mode^2/sqrt((sqrt(pulse.PpeakTheo/medium.P_crit)-0.852)^2-0.0219);
%% Save
save([date,'ION',num2str(beam.n_cycles),'Test.mat']);
                    