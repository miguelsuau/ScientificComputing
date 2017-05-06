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
h = logspace(-2,-1,10);
for i = 1:10
    k(i) = h(i)^2/2;
    mu = 1*k(i)/h(i)^2;
    M = ceil(2/h(i));
    N = ceil(1/k(i));
    U = parabolicSolver(@boundaryFun,h(i),k(i),theta,mu);
    x = linspace(-1,1,M+1);
    t = linspace(0,1,N+1);
    Utrue = zeros(M+1,N);
    for j = 1:N
        Utrue(:,j) = boundaryFun(x,t(j));
    end
    LTE(i) = abs(U(2,2)-Utrue(2,2));
end
%% LTE plots
figure
subplot(1,2,1)
loglog(h,LTE,'o-')
xlabel('h')
ylabel('\tau')
hold on
p1 = loglog(h,h.^2,'--');
legend({'LTE','$\mathcal{O}(h^4)$'},'Interpreter','latex','location','northwest')
subplot(1,2,2)
loglog(k,k,'--')
hold on
loglog(k,LTE,'-o')
xlabel('k')
ylabel('\tau')
legend({'LTE','$\mathcal{O}(k^2)$'},'Interpreter','latex','location','northwest')
print('c2','-dpng')