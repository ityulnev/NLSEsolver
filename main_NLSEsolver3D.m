addpath('Init','help_functions','calc_functions')
clear all
% close all
warning('none')
%% initialize all input Parameters
beam=beam_init;                        % parameters of the beam   #Set(wavelength, energy, duration)  
mesh=mesh_init(beam,1e-3,2);           % create meshgrid          #Set(r,t,f) #In(beam,propagate Distance Lz, Dimension t or t&r)   
medium=medium_init(mesh,beam,'Neon');  % parameters for medium    #Set(refractive index n, n2, beamsize)   
pulse=pulse_init(mesh,beam,medium,0,1);% calculate pulse          #Set(Et(t,r,z=0),Ef,It,If)
%% 2D Propagate, Finite Difference + Split Step
comment='electron drift';
boundcon=["open","open"]; %left and right boundary condition
[Er,Erz,zprop,Qhist,whist,IonizLvl]=do_FourierSplitStep2D(mesh,beam,medium,pulse,boundcon,comment);
%% Self-focusing length
% L_spm=0.367*medium.k0*beam.r_mode^2/sqrt((sqrt(pulse.PpeakTheo/medium.P_crit)-0.852)^2-0.0219);
%% Save
save([date,'ION',num2str(beam.n_cycles),'Test.mat']);

