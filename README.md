# NLSE_2D
Solves the non linear schroedinger like equation for electric field E(t,r) with the low strorage runge kutta method in z.
Transverse component is calculated via finite differences with open boundary conditions.
Valid for multicycle pulses propagating in isotrope media like gases.

Included: Selfphase Modulation + Self focusing, GVD, Divergence, Ionizationloss, Plasma defocusing+ Blue Shift

This version takes as input the envelope with carrier wave iw0t, so the forward propagating complex electric field. With this ionization is calculated dependent for each half cycle of the pulse. This makes the simulation more correct for few-cycle pulses


