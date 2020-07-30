


Iapp = [-100 -70 -40];
it = 0;
r_vec = zeros(3,1);
for l=1:3
	Veff = Iapp(l);
	singlecellmcode_nonoise;
	r_vec(l) = rate*1000;
end


while it<=itmax	%%% bisection for target firing rate

	it = it+1;
	Iapp(2) = (Iapp(3)+Iapp(1))/2;
	Veff = Iapp(2);
	singlecellmcode_nonoise;

	r_vec(2) = rate*1000; % sp/s

	if r_vec(2) > rho0
		r_vec(3) = r_vec(2);
		Iapp(3) = Iapp(2);
	elseif r_vec(2) < rho0
		r_vec(1) = r_vec(2);
		Iapp(1) = Iapp(2);
	else
		break
	end

end

%%% calculate accomodation index
k_disc = 4; % discard this many spikes at the beginning of the trial
accom = 0;

spktimes = spktimes(spktimes>0);
ISI = diff(spktimes);
for j=k_disc:length(ISI)
	accom = accom + (ISI(j)-ISI(j-1))/(ISI(j)+ISI(j-1));
end
accom = accom/(count-k_disc-1);
