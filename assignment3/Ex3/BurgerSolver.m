function U = BurgerSolver(boundaryFun,h,k,epsilon,tmax)
M = ceil(2/h);
N = ceil(tmax/k);
x = linspace(-1,1,M);
t = linspace(0,tmax,N);
U = zeros(M,N);
U(:,1) = boundaryFun(x,0,epsilon);
U(1,:) = boundaryFun(-1,t,epsilon);
U(N,:) = boundaryFun(1,t,epsilon);
idx = 2:(M-1);
for n=2:N
    U(idx,n) = k*epsilon/h^2*(U(idx-1,n-1)-2*U(idx,n-1)+U(idx+1,n-1)) +...
        U(idx,n-1) - k/h*U(idx,n-1).*(U(idx,n-1)-U(idx-1,n-1));
end

end