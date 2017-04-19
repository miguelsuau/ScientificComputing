function [c,info] = bisection(odefun,ts,epsilon,a,b,tol,maxit)
rc = tol+1;
k = 1;

Ua0 = [-1; a];
[ta,Ua] = ode45(odefun,ts,Ua0,[],epsilon);
ra = Ua(end,1)-1.5;

while ((k < maxit) && (abs(rc(end)) > tol))
    c = (a+b)/2;
    Uc0 = [-1; c];
    [tc,Uc] = ode45(odefun,ts,Uc0,[],epsilon);
    rc(k) = Uc(end,1)-1.5;
    if (sign(rc(end)) == sign(ra))
        a = c;
    else
        b = c;
    end
    k = k+1;
end
info.r = abs(rc);
info.iter = k-1;
end