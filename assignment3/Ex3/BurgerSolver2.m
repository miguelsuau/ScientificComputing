function [U,t,x] = BurgerSolver2(h,k,epsilon,tmax)
M = ceil(2/h);
N = ceil(tmax/k);
x = linspace(-1,1,M);
t = linspace(0,tmax,N);
M = length(x);
U = zeros(M,N);
U(:,1) = -sin(pi*x);
U(1,:) = 0;
U(M,:) = 0;
idx = 2:(M-1);
for n=2:N
    U(idx,n) = k*epsilon/h^2*(U(idx-1,n-1)-2*U(idx,n-1)+U(idx+1,n-1)) +...
        U(idx,n-1) - k/h*U(idx,n-1).*(U(idx,n-1)-U(idx-1,n-1));
end
end