close all
clear
%% DIFFUSION PROBLEM
%% Implement the scheme
theta = 0;
h = 0.1;
k = 0.001;
mu = 1*k/h^2;
M = 2/h;
N = 1/k;
U = parabolicSolver(@boundaryFun,h,k,theta,mu);
x = linspace(-1,1,M+1);
t = linspace(0,1,N);
Utrue = zeros(M+1,N);
for j = 1:N
    Utrue(:,j) = boundaryFun(x,t(j));
end
%% Demonstrate convergence
theta = 0;
h = logspace(-2,0,7);
for i = 1:7
    k(i) = h(i)^2/6;
    mu = 1*k(i)/h(i)^2;
    M = ceil(2/h(i));
    N = ceil(1/k(i));
    U = parabolicSolver(@boundaryFun,h(i),k(i),theta,mu);
    x = linspace(-1,1,M+1);
    t = linspace(0,1,N);
    Utrue = zeros(M+1,N);
    for j = 1:N
        Utrue(:,j) = boundaryFun(x,t(j));
    end
    LTE(i) = max(max(abs(U-Utrue)));
end
% LTE plots
figure
subplot(1,2,1)
loglog(h,LTE)
hold on
loglog(h,h.^4,'--')
subplot(1,2,2)
loglog(k,k.^2,'--')
hold on
loglog(k,LTE)