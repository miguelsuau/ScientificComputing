close all
%% 2-POINT BVPs

%% Solve using bvp4c
odefun = @(t,x,epsilon) [x(2); -x(1)*(x(2)-1)/epsilon];
bcfun = @(xa,xb,epsilon) [xa(1)+1; xb(1)-1.5];

options = bvpset('reltol',1e-2,'abstol',1e-2);
options2 = bvpset('reltol',1e-10,'abstol',1e-10);

solinit = bvpinit(linspace(0,1,10),[0 0]);

epsilon = 0.01;
sol = bvp4c(odefun,bcfun,solinit,options,epsilon);
sol2 = bvp4c(odefun,bcfun,solinit,options2,epsilon);
plot(sol2.x,sol2.y(1,:),'linewidth',1.6)
%% Solve using Newton's method
n = length(sol.x);
n = 100;
t = linspace(0,1,n)';
% Aproximate solution 2.105
omega = 0.5*(0-1+1.5+1);
tm = 0.5*(0+1+1-1.5);
Uinit = t-tm+omega*tanh(omega*(t-tm)/(2*0.01));
%plot(t,Uinit,'--','linewidth',1.6)
hold on

    
h = 1/(n+1);
tol = 1e-5;
maxit = 1000;
U = NewtonsMethod(@FunJac,Uinit,h,tol,maxit);
plot(t,U,'linewidth',1.6)

xlabel('t')
ylabel('U')
legend('bvp4c','Newton','location','southeast')
print('sol2','-dpng')
%% Solve using the single shooting method

%% Bisection method
ts = [0 1];
epsilon = 0.01;

a = 1;
Ua0 = [-1; a];
[ta,Ua] = ode45(odefun,ts,Ua0,[],epsilon);

b = 5;
Ub0 = [-1; b];
[tb,Ub] = ode45(odefun,ts,Ub0,[],epsilon);
t1 = cputime;
for i = 1:1
[c,info1] = bisection(odefun,ts,epsilon,a,b,tol,maxit);
end
tbisection = (cputime-t1)/100;

Uc0 = [-1; c];
[tc,Uc] = ode45(odefun,ts,Uc0,[],epsilon);

figure
plot(ta,Ua(:,1),tb,Ub(:,1),'linewidth',1.6)
hold on
plot(tc,Uc(:,1),'linewidth',1.6)

xlabel('t')
ylabel('U')
legend('U_a','U_b','U_c','location','northwest')
print('bisection','-dpng')
%% Sensitivity analysis finite difference approximation
% Forward
epsilon = 0.001;
mu = 1e-7;
U0e = [-1; c+mu];
[te,Ue] = ode45(odefun,ts,U0e,[],epsilon);
duds = (Ue(end,1) - Uc(end,1))/mu;
% Backward
mu = 1e-7;
U0e = [-1; c-mu];
[te,Ue] = ode45(odefun,ts,U0e,[],epsilon);
duds2 = (Uc(end,1) - Ue(end,1))/mu;

% Sensitivity analysis analytical solution
odefun2 = @(t,w,epsilon) [w(2); -w(1)*(w(2)-1)/epsilon; w(4); (-w(2)+1)/epsilon*w(3)-w(1)/epsilon*w(4)];
w0 = [-1; c; 0; 1];
[t,W] = ode45(odefun2,ts,w0,[],epsilon);


%% Plot solution as a function of sigma
sigma = linspace(0.9999,1.01,1000);
Uend = zeros(1000,1);
for i=1:1000   
    U0 = [-1; sigma(i)];    
    [t,U] = ode45(odefun,ts,U0,[],epsilon);
    Uend(i) = U(end,1);
end
r = Uend - 1.5;
subplot(1,2,1)
plot(sigma,r)
xlabel('\sigma')
ylabel('r(\sigma)')
subplot(1,2,2)
plot(sigma,r)
xlabel('\sigma')
ylabel('r(\sigma)')
axis([1 1.00001 -2 0.2]);
print('usigma','-dpng')
%% Newton's method
sigma0 =3;
tol = 1e-2;
maxit = 100;
epsilon = 0.01;
odefun2 = @(t,w,epsilon) [w(2); -w(1)*(w(2)-1)/epsilon; w(4); (-w(2)+1)/epsilon*w(3)-w(1)/epsilon*w(4)];

t1 = cputime;
for i = 1:1000
[sigma,info2] = NewtonShooting(odefun2,ts,epsilon,sigma0,tol,maxit);
end
tnewton = (cputime-t1)/1000;

plot(1:info2.iter,info2.r,'linewidth',1.6)
hold on
plot(1:info1.iter,info1.r,'linewidth',1.6)
axis([1 7 0 0.6])
xlabel('iteration')
ylabel('residual')
legend('Newton','bisection')
print('rate','-dpng')
%% Finite difference approximation/Secant Method
sigma0 = 3;
t1 = cputime;
for i=1:1000
[sigma,info3] = secant(odefun,ts,epsilon,sigma0,tol,maxit);
end
tsecant = (cputime-t1)/1000;