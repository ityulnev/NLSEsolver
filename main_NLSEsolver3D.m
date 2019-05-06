addpath('Init','help_functions','calc_functions')
clear all
close all
%% initialize all input Parameters
beam=beam_init;                      % parameters of the beam   #Set(wavelength, energy, duration)  
mesh=mesh_init(beam,17.1e-3,2);       % create meshgrid          #Set(r,t,f) #In(beam,propagate Distance Lz, Dimension t or t&r)   
medium=medium_init(mesh,beam,'Neon');% parameters for medium    #Set(refractive index n, n2, beamsize)   
pulse=pulse_init(mesh,beam,medium,0,1);% calculate pulse        #Set(Et(t,r,z=0),Ef,It,If)
%% 2D Propagate, Finite Difference + Split Step
comment='Divergence propagation';
boundcon="openCorrected"; %'openCorrected'
for bc=boundcon
[Erf,Erz,zprop,dQhist,whist]=do_FourierSplitStep2D(mesh,beam,medium,pulse,bc,comment);
end

save([date,'Test1.mat']);

