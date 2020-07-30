%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% celltype = 'pyr' % pyr, pv_ or som
Npar = 40;

c = .1;

if length(celltype)==3
	if celltype=='pyr'
		par_pyr;
		R_vec = linspace(2,4.2,Npar); % for pyr
	elseif celltype=='som'
		par_som;	
		R_vec = linspace(1.25,3.1,Npar); % for som
	end
elseif length(celltype)==2
	if celltype=='pv'
		par_pv;
		R_vec = linspace(1.5,3,Npar); % for pv
	end
elseif length(celltype)==4
	if celltype=='test'
		par_test;
		R_vec = linspace(6.5,9,Npar);
	end
elseif length(celltype)==11
	if celltype=='pyr_current'
		par_pyr;
		% RI = 3;
		% RE = 3;
		R_vec = linspace(2,20,Npar);
		% R_vec = linspace(0,2,Npar);
	end
else
	error('set cell type at line 7: pyr, pv or som or pyr_current')
end

T = (1:200);
nT = length(T);

r_vec = zeros(Npar,1);
cov_vec = zeros(Npar,nT);
var_vec = zeros(Npar,nT);
sig_vec = zeros(Npar,1);
tau_vec = zeros(Npar,1);
gain_vec = zeros(Npar,1);

N = 2;

nfreq = 2^12;
l_f = .001/1000;
h_f = 500/1000;
freq = linspace(l_f,h_f,nfreq);

% solve for firing rates via fixed point iteration
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

tau = C/gL;

num_rate_fp_its = 1;

Veff = (VL + tau*RE*g0E*VE + tau*RI*g0I*VI)/(1 + tau*RE*g0E + tau*RI*g0I);
sigma2 = (g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2;

u1 = 1;   
mu = 0;

for np=1:Npar

	RE = R_vec(np);

    %%% diffusion approximation
    if length(celltype)==11
    	if celltype == 'pyr_current'
    		taueff = tau;
    		geff = gL;
    		% mu = R_vec(np);
		    Veff = (VL + tau*RE*g0E*VE + tau*RI*g0I*VI)/(1 + tau*RE*g0E + tau*RI*g0I);
		    sigma2 = (g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2;
    	else
    		error('length(celltype) = 11 but celltype ~= pyr_current');
    	end
    else
	    taueff = tau/(1 + tau*RE*g0E + tau*RI*g0I);
        geff = C/taueff;
        Veff = (VL + tau*RE*g0E*VE + tau*RI*g0I*VI)/(1 + tau*RE*g0E + tau*RI*g0I);
    	sigma2 = (g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2;
	end
    % Veff = (VL + tau*RE*g0E*VE + tau*RI*g0I*VI)/(1 + tau*RE*g0E + tau*RI*g0I);
    % sigma2 = (g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2;


    % Deltaeff = taueff*Delta/tau;

    params(1) = geff;
    params(5) = Veff;
    % params(3) = Deltaeff;

    [P0,p0,~,r0,x0] = theory0(mu,sigma2,params,xi); %initial estimate

	[At,x1,V1] = theory_EIFresp(params,x0,mu,sigma2,xi,u1,r0,P0,p0,freq);
	[f0] = theory_EIFfpt(params,mu,sigma2,freq);
	Ct0 = r0.*(1+2.*real(f0./(1-f0))); %renewal power spectrum

	nT = length(T);
	covT=zeros(nT,1);
	k_T=zeros(nfreq,nT);

	D = (geff/C)*sqrt(2*(C/geff)*sigma2);
	Pss = c*D^2;

	for t=1:nT
	    k_T(:,t)=(sin(pi.*freq.*T(t)).^2)./(pi^2.*freq.^2);
	    cov_vec(np,t)=trapz(freq,Pss.*k_T(:,t).*abs(At).^2);
	    var_vec(np,t)=trapz(freq,(Ct0+Pss.*abs(At).^2).*k_T(:,t));
	end


	r_vec(np) = r0;
	sig_vec(np) = sqrt(sigma2);
	tau_vec(np) = taueff;
	gain_vec(np) = abs(At(1));

end