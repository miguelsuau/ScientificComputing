%% Newton's method Finite Difference Aproximation to the partial derivative
odefun = @(t,x) [x(2); -x(1)*(x(2)-1)/0.01];
tol = 1-3;
maxit = 10000;
sigma = 1.5;
k = 0;
ts = [0; 1];
r = tol+1;
epsilon = 1e-4;
while ((k < maxit) && (abs(r) > tol))
    U0 = [-1; sigma];
    [t1,U1] = ode45(odefun,ts,U0);
    U0e = [-1; sigma+epsilon];
    [t2,U2] = ode45(odefun,ts,U0e);
    duds = (U2(end,1) - U1(end,1))/epsilon;
    r = U1(end,1)-1.5;
    sigma = sigma - 0.001*r/duds;
    k = k+1;
end