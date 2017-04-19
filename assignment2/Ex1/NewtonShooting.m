function [sigma,info] = NewtonShooting(odefun,ts,epsilon,sigma0,tol,maxit)

r = tol+1;
k = 1;
sigma = sigma0;

while ((k < maxit) && (abs(r(end)) > tol))
    w0 = [-1; sigma; 0; 1];
    [t,W] = ode45(odefun,ts,w0,[],epsilon);
    r(k) = W(end,1)-1.5;
    sigma = sigma - r(end)/W(end,3);
    k = k+1;
end

info.r = abs(r);
info.iter = k-1;
end