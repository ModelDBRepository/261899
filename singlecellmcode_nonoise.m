% function [rate] = singlecellmcode(tau, sig, mu, Vth, Vr, tstop, trans, dt);

trans = 1000;
tstop = 1500;
dt = .01;

t = 0;
count = 0;
rate = 0;
tlast = -tstop;
Vs = zeros(tstop/dt,1);
x=zeros(tstop/dt,1);
spktimes = zeros(10000,1);
Vs(1) = Vr;
sqrtdt = sqrt(dt);
g_spike = gL*Delta;

for i=2:tstop/dt
    t=t+dt;
    
    x(i)=x(i-1)+dt*(1./(1+exp(-(Vs(i-1)-Vxh)/Dx))-x(i-1))./tau_x;

    % Vs(i) = Vs(i-1)+(dt/C)*(gL*(VL-Vs(i-1))+gx*x(i)*(Vx-Vs(i-1))+mu+g_spike*exp((Vs(i-1)-VT)/Delta));
    Vs(i) = Vs(i-1)+(dt/taueff)*(Veff-Vs(i-1)+Delta*exp((Vs(i-1)-VT)/Delta))+dt/C*gx*x(i-1)*(Vx-Vs(i-1)); % works with noise



    if t-tlast <= tref
        Vs(i) = Vth-(t-tlast)/tref*(Vth-Vr);
    end
    
    if Vs(i) > Vth
        Vs(i) = Vr;
        tlast = t;
        count = count+1;
        spktimes(count)=t;

        if t>trans
            rate = rate+1;
        end
    end

end

rate = rate./(tstop-trans);
spktimes = spktimes(spktimes>0);