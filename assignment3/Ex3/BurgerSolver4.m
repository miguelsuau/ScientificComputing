function U = BurgerSolver4(boundaryFun,h,k,epsilon,tmax)
M = ceil(2/h);
N = ceil(tmax/k);
x = linspace(-1,1,M+1);
t = linspace(0,tmax,N+1);
U = zeros(M+1,N+1);
U(:,1) = boundaryFun(x,0,epsilon)';
U(1,:) = boundaryFun(-1,t,epsilon);
U(M,:) = boundaryFun(1,t,epsilon);
idx = 2:(M-1);
for n=2:N+1
    fm = 0.5*(U(2:M-1,n-1).^2 + U(1:M-2,n-1).^2);
    fp = 0.5*(U(2:M-1,n-1).^2 + U(3:M,n-1).^2);
    U(idx,n) = k*epsilon/h^2*(U(idx-1,n-1)-2*U(idx,n-1)+U(idx+1,n-1)) +...
        U(idx,n-1) - k/(2*h)*(fp-fm);
end
end