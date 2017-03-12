addpath('../')

t0 = 0;
tf = 50;
x0 = [2; 0];
mu = 3; % and then 100
h0 = 0.01;
absTol = 1e-6;
relTol = 1e-6;


% Adaptive ESDIRK23
[Tout,Xout,info,stats] = ESDIRK23_Adaptive(@VanderPolFun,@VanderPolJac,t0,tf,x0,h0,absTol,relTol,mu);
t = Tout;
x = Xout;

% plot results
figure
subplot(4,2,1)
plot(t,x(:,1),'.-b'), title('Van der Pol solution - ESDIRK_{23}^{Adaptive} \mu = 3'), xlim([t0 tf])
subplot(4,2,2)
plot(t,x(:,2),'.-r'), xlabel('Time'), xlim([t0 tf])

% ESDIRK23
[Tout,Xout,info,stats] = ESDIRK23(@VanderPolFun,@VanderPolJac,t0,tf,x0,h0,absTol,relTol,mu);
t2 = Tout;
x2 = Xout;

subplot(4,2,3)
plot(t2,x2(:,1),'.-b'), title('Van der Pol solution - ESDIRK_{23}^{Fixed} \mu = 3'), xlim([t0 tf])
subplot(4,2,4)
plot(t2,x2(:,2),'.-r'), xlabel('Time'), xlim([t0 tf])

% change mu to 100
mu = 100;
h0 = 0.01;
tf = 1000;

[Tout,Xout,info,stats] = ESDIRK23_Adaptive(@VanderPolFun,@VanderPolJac,t0,tf,x0,h0,absTol,relTol,mu);
t = Tout;
x = Xout;
subplot(4,2,5)
plot(t,x(:,1),'.-b'), title('Van der Pol solution - ESDIRK_{23}^{Adaptive} \mu = 100'), xlim([t0 tf])
subplot(4,2,6)
plot(t,x(:,2),'.-r'), xlabel('Time'), xlim([t0 tf])

[Tout,Xout,info,stats] = ESDIRK23(@VanderPolFun,@VanderPolJac,t0,tf,x0,h0,absTol,relTol,mu);
t2 = Tout;
x2 = Xout;
subplot(4,2,7)
plot(t2,x2(:,1),'.-b'), title('Van der Pol solution - ESDIRK_{23}^{Fixed} \mu = 100'), xlim([t0 tf])
subplot(4,2,8)
plot(t2,x2(:,2),'.-r'), xlabel('Time'), xlim([t0 tf])