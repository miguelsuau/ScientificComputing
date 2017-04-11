close all
%% 2-POINT BVPs

%% Solve using bvp4c
odefun = @(t,x) [x(2); -x(1)*(x(2)-1)/0.01];
bcfun = @(xa,xb) [xa(1)+1; xb(1)-1.5];
options = bvpset('reltol',1e-2,'abstol',1e-2);
solinit = bvpinit(linspace(0,1,10),[0 0]);
sol = bvp4c(odefun,bcfun,solinit);%,options);

plot(sol.x,sol.y(1,:))
%% Solve using Newton's method
t = linspace(0,1,1000)';
% Aproximate solution 2.105
omega = 0.5*(0-1+1.5+1);
tm = 0.5*(0+1+1-1.5);
Uinit = t-tm+omega*tanh(omega*(t-tm)/(2*0.01));
plot(t,Uinit)
hold on
    
h = 1/999;
tol = 1e-5;
maxit = 1000;
U = NewtonsMethod(@FunJac,Uinit,h,tol,maxit);
plot(t,U)

%% Solve using the single shooting method

%% Bisection method
k = 0;
rc = tol+1;
ts = [0 1];

a = 1;
Ua0 = [-1; a];
[ta,Ua] = ode45(odefun,ts,Ua0);
ra = Ua(end,1)-1.5;

b = 10;
Ub0 = [-1; b];
[tb,Ub] = ode45(odefun,ts,Ub0);

figure
plot(ta,Ua(:,1),tb,Ub(:,1))
hold on
while ((k < maxit) && (abs(rc) > tol))
    c = (a+b)/2;
    Uc0 = [-1; c];
    [tc,Uc] = ode45(odefun,ts,Uc0);
    rc = Uc(end,1)-1.5;
    if (sign(rc) == sign(ra))
        a = c;
    else
        b = c;
    end
    k = k+1;
end
plot(tc,Uc(:,1))

% Sensitivity analysis finite difference approximation
epsilon = 1e-7;
U0e = [-1; c+epsilon];
[te,Ue] = ode45(odefun,ts,U0e);
duds = (Ue(end,1) - Uc(end,1))/epsilon;

% Sensitivity analysis analytical solution
odefun2 = @(t,w) [w(2); -w(1)*(w(2)-1)/0.01; w(4); (-w(2)+1)/0.01*w(3)-w(1)/0.01*w(4)];
w0 = [-1; c; 0; 1];
[t,W] = ode45(odefun2,ts,w0);


%% Newton's method
sigma = 5;
k = 0;
tol = 1e-2;
r = tol+1;
odefun2 = @(t,w) [w(2); -w(1)*(w(2)-1)/0.01; w(4); (-w(2)+1)/0.01*w(3)-w(1)/0.01*w(4)];
while ((k < maxit) && (abs(r) > tol))
    w0 = [-1; sigma; 0; 1];
    [t,W] = ode45(odefun2,ts,w0);
    r = W(end,1)-1.5;
    sigma = sigma - r/W(end,3);
    k = k+1;
end


%% Finite difference approximation/Secant Method
sigma = 1;
k = 0;
r = tol+1;
epsilon = 1e-4;
while ((k < maxit) && (abs(r) > tol))
    U0 = [-1; sigma];
    [t1,U1] = ode45(odefun,ts,U0);
    U0e = [-1; sigma+epsilon];
    [t2,U2] = ode45(odefun,ts,U0e);
    duds = (U2(end,1) - U1(end,1))/epsilon;
    r = U1(end,1)-1.5;
    sigma = sigma - r/duds;
    k = k+1;
end