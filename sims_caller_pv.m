%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% call sims, looping over firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mex -g CFLAGS="\$CFLAGS -std=c99" spk_cond.c;

celltype = 'pv' % pyr, pv or som

if length(celltype)==3 && celltype=='pyr'
	par_pyr;
	R_vec = linspace(4,5.5,Npar); % for pyr
elseif length(celltype)==2 && celltype == 'pv'
	par_PV;
	R_vec = linspace(1,2,Npar); % for pv
elseif length(celltype)==3 && celltype == 'som'
	par_SOM;	
	R_vec = linspace(5.5,9,Npar); % for som
else
	error('set cell type at line 7: pyr, pv or som')
end

Npar = 50;

T = 2:2:200;
nT = length(T);

r_vec = zeros(Npar,1);
cov_vec = zeros(Npar,nT);
var_vec = zero(Npar,nT);

if matlabpool('size')==0
	matlabpool open 10;
end

for np = 1:Npar

	RE = R_vec(np);
	taueff = tau/(1 + tau*RE*g0E + tau*RI*g0I);
	Veff = (VL + tau*RE*g0E*VE + tau*RI*g0I*VI)/(1 + tau*RE*g0E + tau*RI*g0I);
	sig2eff = (g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2;
	sigeff = sqrt(sig2eff);

	sims_pair_cond;

	r_vec(np) = rateE;
	cov_vec(np,:) = cov_calc;
	var_vec(np,:) = var_calc;

end

if matlabpool('size')~=0
	matlabpool close
end

save(strcat('sims_',celltype,'_c=',num2str(c),'_',datestr(now,1)));

exit