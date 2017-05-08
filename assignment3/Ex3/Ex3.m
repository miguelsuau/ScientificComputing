%% Non-linear advection-diffusion equation
h = 0.005;
k = h^2;
epsilon = 0.5;
tmax = 1;
U1 = BurgerSolver(@boundaryFun,h,k,epsilon,tmax);
%% true solution
figure
x = linspace(-1,1,100);
t = linspace(0,1,100);
[X,T] = ndgrid(x,t);
U2 = boundaryFun(X,T,epsilon);
mesh(X,T,U2)

%% Demonstrate convergence

h = [1/50 1/100 1/150 1/200];
k = h.^2;
for i = 1:4
    theta(i) = 0.5 + h(i)^2/(12*k(i));
    mu = 1*k(i)/h(i)^2;
    M = ceil(2/h(i));
    N = ceil(1/k(i));
    U =BurgerSolver4(@boundaryFun,h(i),k(i),epsilon,tmax);
    x = linspace(-1,1,M+1);
    t = linspace(0,1,N+1);
    [X,T] = ndgrid(x,t);
    Utrue = boundaryFun(X,T,epsilon);
    LTE(i) = max(norm(U(2,2)-Utrue(2,2),2));
end

%% LTE plots
figure
subplot(1,2,1)
loglog(h,LTE,'o-')
xlabel('h')
ylabel('\tau')
hold on
p1 = loglog(h,h.^3,'--');
legend({'LTE','$\mathcal{O}(h^3)$'},'Interpreter','latex','location','southeast')
subplot(1,2,2)
loglog(k,LTE,'-o')
hold on
loglog(k,k.^2*10,'--')
xlabel('k')
ylabel('\tau')
legend({'LTE','$\mathcal{O}(k^2)$'},'Interpreter','latex','location','southeast')
print('c3','-dpng')
%% 3.3

h = 0.005;
k = 1/2*h^2;
tmax = 1.6037/pi;
xspan = [-1 1];
epsilon = 0.01/pi;
%% FTCS and upwind
[U1,t1,x1] = BurgerSolver2(h,k,epsilon,tmax);
plot(x1,U1(:,end))
% find closest point to 0
[~,idx] = min(abs(x2));
dx1 = (U(idx+1,end)-U(idx,end))/h;
%% Trapezoidal
[U2,t2,x2] = BurgerSolver5(h,k,epsilon,tmax);
[T,X] = meshgrid(t,x);
plot(x2,U2(:,end))
% find closest point to 0
[~,idx] = min(abs(x2));
dx2 = (U(idx+1,end)-U(idx,end))/h;
%% Nonuniform grid
% Divide the problem in sections
Xspan = [-1 -0.09];
h = 0.1;
k = h^2;
[U1,t1,x1] = BurgerSolver5(h,k,epsilon,tmax,xspan);
xpan = [-0.1 0.1];
h = 0.01;
k = h^2;
[U2,t2,x2] = BurgerSolver5(h,k,epsilon,tmax,xspan);
xspan = [0.11 1];
h = 0.1;
k = h^2;
[U3,t3,x3] = BurgerSolver5(h,k,epsilon,tmax,xspan);
% join solution