close all
clear
%% DIFFUSION PROBLEM
%% Implement the scheme
theta = 0.5;
h = 0.01;    
k = 0.01;
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
subplot(2,2,1)
plot(x,U(:,1),'linewidth',1.6)
title('t=0')
subplot(2,2,2)
plot(x,U(:,5),'linewidth',1.6)
title('t=0.04')
subplot(2,2,3)
plot(x,U(:,20),'linewidth',1.6)
title('t=0.20')
subplot(2,2,4)
plot(x,U(:,50),'linewidth',1.6)
title('t=0.5')
print('c3','-dpng')
%% Demonstrate convergence

h = logspace(-2,-1,10);
k = 7/6*h.^2;
for i = 1:10
    theta(i) = 0.5 + h(i)^2/(12*k(i));
    mu = 1*k(i)/h(i)^2;
    M = ceil(2/h(i));
    N = ceil(1/k(i));
    U = parabolicSolver(@boundaryFun,h(i),k(i),theta(i),mu);
    x = linspace(-1,1,M+1);
    t = linspace(0,1,N+1);
    Utrue = zeros(M+1,N);
    for j = 1:2
        Utrue(:,j) = boundaryFun(x,t(j));
    end
    LTE(i) = max(norm(U(:,2)-Utrue(:,2),2));
end

%% LTE plots
figure
subplot(1,2,1)
loglog(h,LTE,'o-')
xlabel('h')
ylabel('\tau')
hold on
p1 = loglog(h,h.^2,'--');
legend({'LTE','$\mathcal{O}(h^2)$'},'Interpreter','latex','location','northwest')
subplot(1,2,2)
loglog(k,LTE,'-o')
hold on
loglog(k,k.^2,'--')
xlabel('k')
ylabel('\tau')
legend({'LTE','$\mathcal{O}(k^2)$'},'Interpreter','latex','location','southeast')
print('c3','-dpng')