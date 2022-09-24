% This program plots the dispersion relation h(k^2) in the Schnakenberg
% reaction-diffusion model. 
% Mathematical Biology modeling project, 2018 Michaelmas Term, Oxford
% University. 

Du = 1; Dv = 40; 
uu = 0.9; vv = 1;
fu = @(u,v) - 1 + 2*u*v;
fv = @(u,v) u^2; 
gu = @(u,v) -2*u*v;
gv = @(u,v) -u^2;

h = @(ksq,uu,vv) Du*Dv*ksq.^2 - (Dv*fu(uu,vv) + Du*gv(uu,vv)).*ksq + (fu(uu,vv)*gv(uu,vv) - gu(uu,vv)*fv(uu,vv))
Ksq = linspace(0,1,50);
plot(Ksq,zeros(1,length(Ksq)),'linewidth',2); hold on
plot(Ksq,h(Ksq,uu,vv),'linewidth',3); xlabel('k^2'); ylabel('h(k^2)');
set(gca,'fontsize',20)
grid on