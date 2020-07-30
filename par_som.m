%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cellular and input parameters - putative SOM+ cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NE = 2;

%%% cellular parameters
gL = .04; %mS/sq cm
% VT = -46;
VT = -52;
Vr = -65;
% Vr = -;
Vth = 0;
VL = -65; % mV
% VL = -50;	
C = 1; %uF/sq cm^2
tref = 4;
tau = C/gL;
% tau = 25;
Delta = 1.4;

Vlb = -100;
dV = .1;
V = Vlb+dV:dV:Vth;

gx = .1; %mS/sq cm
Dx = 8; %slope factor
Vx = -85; %reversal potential of KCNQ
Vxh = -40; %half-activation of KCNQ, -45 / -37
xi = 1./(1+exp(-(V-Vxh)/Dx)); %steady-state activation x_infinity
tau_x = 200; %ms, only measured for control
% tau_x = 100;

%%% input parameters
% g0E = .02*gL; % strength of excitatory kicks
% g0I = .06*gL; % strength of inhibitory kicks
% g0E = .01;
% g0I = .02;
% VE = 0; % excitatory reversal potential
% VI = -80; % inhibitory reversal potential
% RI = 6; %1/ms, fixed
% RE = 9; 

g0E = .06;
g0I = .06;
VE = 0; % excitatory reversal potential
VI = -80; % inhibitory reversal potential
RI = 10; %1/ms, fixed
RE = 3; 

%%% diffusion approximation
taueff = tau/(1 + tau*RE*g0E + tau*RI*g0I);
Veff = (VL + tau*RE*g0E*VE + tau*RI*g0I*VI)/(1 + tau*RE*g0E + tau*RI*g0I);
sig2eff = (g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2;
sigeff = sqrt(sig2eff);