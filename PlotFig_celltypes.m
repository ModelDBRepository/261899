%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot figure: theory and sims comparing 3 cell types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% compile C codes
mex thin_x.c;
mex theory_EIFresp.c;
mex theory_EIFfpt.c;
mex spk_cond.c;

figure;

%%% voltage traces
clear all;
par_pyr;
RI = 10;
RE = 3;
taueff = tau/(1 + tau*RE*g0E + tau*RI*g0I);
Veff = (VL + tau*RE*g0E*VE + tau*RI*g0I*VI)/(1 + tau*RE*g0E + tau*RI*g0I);
sigma2 = (g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2;
sigeff = sqrt(sigma2);
geff = C/taueff;
mu = 0;

singlecellmcode;
tplot = (dt:dt:tstop)-trans;
subplot(6,3,1); 
plot(dt:dt:tstop,Vs,'k');
axis([0 500 -90 0]); box off; set(gca,'XTick',[]);

clear all;
par_pv;
RI = 10;
RE = 3;
taueff = tau/(1 + tau*RE*g0E + tau*RI*g0I);
Veff = (VL + tau*RE*g0E*VE + tau*RI*g0I*VI)/(1 + tau*RE*g0E + tau*RI*g0I);
sigma2 = (g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2;
sigeff = sqrt(sigma2);
geff = C/taueff;
mu = 0;

singlecellmcode;
tplot = (dt:dt:tstop)-trans;
subplot(6,3,4); plot(tplot,Vs,'b');
axis([0 500 -90 0]); box off; set(gca,'XTick',[]);

clear all;
par_som;
RI = 10;
RE = 3;
taueff = tau/(1 + tau*RE*g0E + tau*RI*g0I);
Veff = (VL + tau*RE*g0E*VE + tau*RI*g0I*VI)/(1 + tau*RE*g0E + tau*RI*g0I);
sigma2 = (g0E^2)*RE*(VE-Veff)^2 + (g0I^2)*RI*(VI-Veff)^2;
sigeff = sqrt(sigma2);
geff = C/taueff;
mu = 0;

singlecellmcode;
tplot = (dt:dt:tstop)-trans;
subplot(6,3,7); plot(tplot,Vs,'r');
axis([0 500 -90 0]); box off;


%%% Bar graphs: VL, VT, passive membrane time constant, accomodation indices, rheobase
clear all;
par_pyr;
subplot(2,15,6);
bar(1,VL,'facecolor','k','edgecolor','none');
hold on; box off;
axis([0 4 -85 -55]);
title('V_L');

subplot(2,15,7);
bar(1,VT,'facecolor','k','edgecolor','none');
hold on; box off;
axis([0 4 -55 -40]);
title('V_T');

subplot(2,15,8);
bar(1,tau,'facecolor','k','edgecolor','none');
hold on; box off;
axis([0 4 0 45]);
title('\tau');

subplot(2,15,9);
rho0 = 15;
itmax = 20;
accomodation_index;
bar(1,accom,'facecolor','k','edgecolor','none');
hold on; box off;
axis([0 4 0 .05]);
title('Accom. Index');

subplot(2,15,10);
clear all; par_pyr;
theory_fIcurve_rheobase;
Ibase = Ivec(min(find(r_vec*1000 > .1)));
bar(1,Ibase,'facecolor','k','edgecolor','none');
hold on; box off;
axis([0 4 0 1]);
title('Rheobase');

clear all;
par_pv;
subplot(2,15,6);
bar(2,VL,'facecolor','b','edgecolor','none');

subplot(2,15,7);
bar(2,VT,'facecolor','b','edgecolor','none');

subplot(2,15,8);
bar(2,tau,'facecolor','b','edgecolor','none');

subplot(2,15,9);
rho0 = 50;
itmax = 20;
accomodation_index;
bar(2,accom,'facecolor','b','edgecolor','none');

subplot(2,15,10);
clear all; par_pv;
theory_fIcurve_rheobase;
Ibase = Ivec(min(find(r_vec*1000 > .1)));
bar(2,Ibase,'facecolor','b','edgecolor','none');

clear all;
par_som;
subplot(2,15,6);
bar(3,VL,'facecolor','r','edgecolor','none');

subplot(2,15,7);
bar(3,VT,'facecolor','r','edgecolor','none');

subplot(2,15,8);
bar(3,tau,'facecolor','r','edgecolor','none');

subplot(2,15,9);
rho0 = 20;
itmax = 20;
accomodation_index;
bar(3,accom,'facecolor','r','edgecolor','none');

subplot(2,15,10);
clear all; par_som;
theory_fIcurve_rheobase;
Ibase = Ivec(min(find(r_vec*1000 > .1)));
bar(3,Ibase,'facecolor','r','edgecolor','none');


%%% f-I curves
clear all;
par_pyr;
theory_fIcurve;
Ibase = Ivec(min(find(r_vec*1000 > .1)));
subplot(2,3,3); plot(Ivec-Ibase,r_vec*1000,'k','linewidth',2);
% axis([-.5 2.5 0 60]); hold on; box off;
axis([0 18 0 150]); hold on; box off;

clear all;
par_pv;
theory_fIcurve;
Ibase = Ivec(min(find(r_vec*1000 > .1)));
subplot(2,3,3); plot(Ivec-Ibase,r_vec*1000,'b','linewidth',2);

clear all;
par_som;
theory_fIcurve;
Ibase = Ivec(min(find(r_vec*1000 > .1)));
subplot(2,3,3); plot(Ivec-Ibase,r_vec*1000,'r','linewidth',2);
legend('Pyr','PV+','SOM+');
xlabel('Current from rheobase (\muA/cm^2)'); ylabel('Firing Rate (sp/s)');


%%% f-I curve gain, effective noise intensity and spike count covariance as functions of rate
%%% first, theory
clear all;
celltype  = 'pyr';
theory_caller;

subplot(2,3,4);
plot(r_vec*1000,gain_vec*1000,'k','linewidth',2);
hold on; box off;
axis([0 60 0 15]); xlabel('Firing Rate (sp/s)'); ylabel('Response Gain ((sp/s) / (\muA/cm^2))');

subplot(2,3,5);
plot(r_vec*1000,sig_vec,'k','linewidth',2);
hold on; box off;
axis([0 60 4 9]); xlabel('Firing Rate (sp/s)'); ylabel('Noise Intensity (mV)');


subplot(2,3,6);
plot(r_vec*1000,cov_vec(:,end),'k','linewidth',2);
hold on; box off;
axis([0 60 0 .21]); xlabel('Firing Rate (sp/s)'); ylabel('Spike Count Cov., T=200 ms (sp/s)^2');


clear all;
celltype  = 'pv';
theory_caller;

subplot(2,3,4);
plot(r_vec*1000,gain_vec*1000,'b','linewidth',2);

subplot(2,3,5);
plot(r_vec*1000,sig_vec,'b','linewidth',2);

subplot(2,3,6);
plot(r_vec*1000,cov_vec(:,end),'b','linewidth',2);

clear all;
celltype  = 'som';
theory_caller;

subplot(2,3,4);
plot(r_vec*1000,gain_vec*1000,'r','linewidth',2);

subplot(2,3,5);
plot(r_vec*1000,sig_vec,'r','linewidth',2);

subplot(2,3,6);
plot(r_vec*1000,cov_vec(:,end),'r','linewidth',2);

%%% second, sims
%%% if have sims saved - do this
% clear all;
% load sims_pyr_c=0.1_06-Jul-2015.mat;
% subplot(2,3,6);
% plot(r_vec,cov_vec(:,end),'ko');

% load sims_pv_c=0.1_06-Jul-2015.mat;
% subplot(2,3,6);
% plot(r_vec,cov_vec(:,end),'bo');

% load sims_som_c=0.1_06-Jul-2015.mat;
% subplot(2,3,6);
% plot(r_vec,cov_vec(:,end),'ro');

%%% if running sims - do this (will take some time obviously, depending on computer)
% clear all;
% sims_caller_pyr;
% subplot(2,3,6); hold on;
% plot(r_vec,cov_vec(:,end),'ko');

% clear all;
% sims_caller_pv;
% subplot(2,3,6); hold on;
% plot(r_vec,cov_vec(:,end),'bo');

% clear all;
% sims_caller_som;
% subplot(2,3,6); hold on;
% plot(r_vec,cov_vec(:,end),'ro');

pos=[0 0 17 10];
set(gcf,'units','centimeters','paperunits','centimeters','PaperPositionMode','Auto','pos',pos,'papersize',[pos(3) pos(4)]);

