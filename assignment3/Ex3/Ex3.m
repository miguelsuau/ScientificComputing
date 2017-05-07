%% Non-linear advection-diffusion equation
h = 0.005;
k = h^2;
epsilon = 0.5;
tmax = 1;
U1 = BurgerSolver4(@boundaryFun,h,k,epsilon,tmax);
mesh(U1);
%% true solution
figure
x = linspace(-1,1,100);
t = linspace(0,1,100);
[X,T] = ndgrid(x,t);
U2 = boundaryFun(X,T,epsilon);
mesh(X,T,U2)
%% 3.3

h = 0.003;
k = h^2;
tmax = 1.6037/pi;
xspan = [-1 1];
epsilon = 0.01/pi;
%% FTCS and upwind
[U1,t1,x1] = BurgerSolver2(h,k,epsilon,tmax);
plot(x1,U1(:,end))
% find closest point to 0
[~,idx] = min(abs(x2));
dx1 = (U(idx+1,end)-U(idx,end))/h;
[T,X] = meshgrid(t1,x1);
mesh(X,T,U1)
%% Trapezoidal
[U2,t2,x2] = BurgerSolver5(h,k,epsilon,tmax);
[T,X] = meshgrid(t,x);
figure
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