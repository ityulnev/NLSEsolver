# NLSE_2D
Solves the non linear schroedinger like equation for electric field E(t,r) with the low strorage runge kutta method in z.
Transverse component is calculated via finite differences with open boundary conditions.
Valid for multicycle pulses propagating in isotrope media like gases.

Included: Selfphase Modulation + Self focusing, GVD, Divergence, Ionizationloss, Plasma defocusing+ Blue Shift

This version 'No_Carrier' takes input electric field without the carrier i*w0*t oscillation to reduce the the fineness of mesh in time. This reduces the arrays sizes in f and t drastically and improves performance. 




