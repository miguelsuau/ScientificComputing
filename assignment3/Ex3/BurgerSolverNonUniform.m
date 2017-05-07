function [U,t,x] = BurgerSolverNonUniform(xNumPts, k,epsilon,tmax)
% Not sure how to deal with nonuniform grid so the implementation is
% adjusted to a bit so it doesn't blow up when deviding by zeros.
%
M = xNumPts; % 667 for example
N = ceil(tmax/k);
x = nonuniformgrid(M);
t = linspace(0,tmax,N);
M = length(x);
U = zeros(M,N);
U(:,1) = -sin(pi*x);
U(1,:) = 0;
U(M,:) = 0;
idx = 2:(M-1);
h = 0.01; % 0.01 for example

diffx = diff(x);
diffx(diffx < 1e-2) = 0.009; % replace small differences with a value for numerical stability

for n=2:N
    fm = 0.5*(U(2:M-1,n-1).^2 + U(1:M-2,n-1).^2);
    fp = 0.5*(U(2:M-1,n-1).^2 + U(3:M,n-1).^2);
    %U(idx,n) = k*epsilon/h^2*(U(idx-1,n-1)-2*U(idx,n-1)+U(idx+1,n-1)) +...
    %    U(idx,n-1) - k/(2*h)*(fp-fm);
    M1 = k*epsilon*diag(diffx.^-2); %k*epsilon/h^2;
    M2 = 0.5*k*diag(diffx.^-1);%k/(2*h);
    U(idx,n) = M1(1:end-1,1:end-1)*(U(idx-1,n-1)-2*U(idx,n-1)+U(idx+1,n-1)) +...
        U(idx,n-1) - M2(1:end-1,1:end-1)*(fp-fm);
end
end