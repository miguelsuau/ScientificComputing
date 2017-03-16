close all; clear all;
addpath('../')
addpath('../Ex1')

f = @(x) [  x(1)+x(2)+x(3)-1;
            1/4*x(2)+x(3)-1/2;
            1/16*x(2)+x(3)-1/3;
            1/4*x(3)*x(6)-1/6;
            x(4)-1/4;
            x(5)+x(6)-1;
         ];
x = fsolve(f,zeros(6,1));
b1 = x(1); b2 = x(2); b3 = x(3); a21 = x(4); a31 = x(5); a32 = x(6);

butcher.AT = [0 0 0; a21 0 0; a31 a32 0]';
butcher.b  = [b1 b2 b3]';
butcher.c  = [0 1/3 1]';
butcher.d  = [1/8 1/2 3/8]' - [b1 b2 b3]'; % d  = bhat-b;
butcher.bhat = [1/8 1/2 3/8]';
butcher.stages = 3;

t0 = 0;
tf = 20;
x0 = [2; 0];
mu = 3; % and then 100
h0 = 0.001;
absTol = 1e-6;
relTol = 1e-6;


% ERK
n1 = ceil((tf - t0)/h0 + 1);
[Tout,Xout,~,info1] = ExplicitRungeKutta(@VanderPolfunjac,[t0;tf],x0,n1,butcher,mu);
t = Tout;
x = Xout;

% ESDIRK23
[Tout,Xout,info2,stats2] = ESDIRK23(@VanderPolFun,@VanderPolJac,t0,tf,x0,h0,absTol,relTol,mu);
t2 = Tout;
x2 = Xout;

subplot(2,2,1)
plot(t2,x2(:,1),'.b'), title('Van der Pol solution \mu = 3');
xlim([t0 tf]), hold on, xlabel('Time'), ylabel('x_1')
plot(t,x(:,1),'.g')
legend('ESDIRK','ERK')

subplot(2,2,2)
plot(t2,x2(:,2),'.b'), xlabel('Time'), ylabel('x_2'), xlim([t0 tf]), hold on
plot(t,x(:,2),'.g'), title('Van der Pol solution \mu = 3')
legend('ESDIRK','ERK')

% change mu to 100 (stiff)
mu = 100;
tf = 200;

% ERK
n2 = ceil((tf - t0)/h0 + 1);
[Tout,Xout,~,info3] = ExplicitRungeKutta(@VanderPolfunjac,[t0;tf],x0,n2,butcher,mu);
t = Tout;
x = Xout;

% ESDIRK23
[Tout,Xout,info4,stats4] = ESDIRK23(@VanderPolFun,@VanderPolJac,t0,tf,x0,h0,absTol,relTol,mu);
t2 = Tout;
x2 = Xout;

subplot(2,2,3)
plot(t2,x2(:,1),'.b'), title('Van der Pol solution \mu = 100');
xlim([t0 tf]), hold on, xlabel('Time'), ylabel('x_1')
plot(t,x(:,1),'.g')
legend('ESDIRK','ERK')

subplot(2,2,4)
plot(t2,x2(:,2),'.b'), xlabel('Time'), ylabel('x_2'), xlim([t0 tf]), hold on
plot(t,x(:,2),'.g'), title('Van der Pol solution - ERK \mu = 100')
legend('ESDIRK','ERK')