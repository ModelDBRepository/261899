%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate f-I curve with fixed noise intensity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% weak noisy current
RI = 2;
RE = 2;
taueff = tau/(1 + tau*RE*g0E + tau*RI*g0I);
Veff = (VL + tau*RE*g0E*VE + tau*RI*g0I*VI)/(1 + tau*RE*g0E + tau*RI*g0I);
sigma2 = (g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2;
sigeff = sqrt(sigma2);
geff = C/taueff;


Ivec = -1:.05:4;
% Ivec = -10:.1:25;
Np = length(Ivec);


params = zeros(13,1);
params(1) = gL;
params(2) = C;
params(3) = Delta;
params(4) = VT;
params(5) = VL;
params(6) = Vth;
params(7) = Vlb;
params(8) = dV;
params(9) = Vr;
params(10) = tref;
params(11) = tau_x;
params(12) = Vx;
params(13) = gx;

r_vec = zeros(Np,1);

for np = 1:Np

	% Veff = Ivec(np);
	% params(5) = Veff;

    [P0,p0,~,r0,x0] = theory0(Ivec(np),sigma2,params,xi); %initial estimate
    r_vec(np) = r0;

end