addpath('../')

t0 = 0;
tf = 50;
x0 = [2; 0];
n = 10000;
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
subplot(2,2,1)
plot(t,x(:,1),'.-b'), title('Van der Pol solution')
subplot(2,2,2)
plot(t,x(:,2),'.-r'), xlabel('Time')

% ESDIRK23
[Tout,Xout,info,stats] = ESDIRK23(@VanderPolFun,@VanderPolJac,t0,tf,x0,h0,absTol,relTol,mu);
t2 = Tout;
x2 = Xout;

subplot(2,2,3)
plot(t2,x2(:,1),'.-b'), title('Van der Pol solution')
subplot(2,2,4)
plot(t2,x2(:,2),'.-r'), xlabel('Time')