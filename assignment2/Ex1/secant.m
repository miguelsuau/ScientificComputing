function [sigma,info] = secant(odefun,ts,epsilon,sigma0,tol,maxit)

k = 1;
r = tol+1;
mu = 1e-4;
sigma = sigma0;

while ((k < maxit) && (abs(r(end)) > tol))
    U0 = [-1; sigma];
    [t1,U1] = ode45(odefun,ts,U0,[],epsilon);
    U0e = [-1; sigma+mu];
    [t2,U2] = ode45(odefun,ts,U0e,[],epsilon);
    duds = (U2(end,1) - U1(end,1))/mu;
    r(k) = U1(end,1)-1.5;
    sigma = sigma - r(end)/duds;
    k = k+1;
end
info.r = abs(r);
info.iter = k-1;
end