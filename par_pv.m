%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cellular and input parameters - putative PV+ cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NE = 2;

%%% cellular parameters
gL = .0667; %mS/sq cm
% VT = -42;
VT = -50;
Vr = -80;
Vth = 0;
% Vth = VT;
% VL = -65; % mV
VL = -80;
C = 1; %uF/sq cm^2
tref = 1;
tau = C/gL;
% tau = 15;
% tau = 25;
% Delta = 1.4;
Delta = .2;

Vlb = -100;
dV = .1;
V = Vlb+dV:dV:Vth;

gx = 0; %mS/sq cm
Dx = 8; %slope factor
Vx = -85; %reversal potential of KCNQ
Vxh = -40; %half-activation of KCNQ, -45 / -37
xi = 1./(1+exp(-(V-Vxh)/Dx)); %steady-state activation x_infinity
tau_x = 200; %ms, only measured for control

%%% input parameters
% g0E = .02*gL; % strength of excitatory kicks
% g0I = .06*gL; % strength of inhibitory kicks
% g0E = .01;
% g0I = .02;
% VE = 0; % excitatory reversal potential
% VI = -80; % inhibitory reversal potential
% RI = 3; %1/ms, fixed
% RE = 8; 

g0E = .06;
g0I = .06;
VE = 0; % excitatory reversal potential
VI = -80; % inhibitory reversal potential
RI = 10; %1/ms, fixed
RE = 8; 


%%% diffusion approximation for effective time constant
taueff = tau/(1 + tau*RE*g0E + tau*RI*g0I);
Veff = (VL + tau*RE*g0E*VE + tau*RI*g0I*VI)/(1 + tau*RE*g0E + tau*RI*g0I);
sig2eff = ((g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2);

% sig2eff = ((g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2);

sigeff = sqrt(sig2eff);
Deltaeff = taueff*Delta/tau;
geff = C/taueff;

bE = 1-exp(-g0E); 
bI = 1-exp(-g0I); 
% bE = bE+bE^2;
% bI = bI+bI^2;

% taueff = tau/(1 + tau*RE*bE + tau*RI*bI);
% Veff = (VL + tau*RE*bE*VE + tau*RI*bI*VI)/(1 + tau*RE*bE + tau*RI*bI);
% sig2eff = (bE^2)*RE*(VE-Veff)^2 + (bI^2)*RI*(VI-Veff)^2;
% sigeff = sqrt(sig2eff);
% Deltaeff = taueff*Delta/tau;
