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
h = 0.005;
k = h^2;
tmax = 1.6037/pi;
epsilon = 0.01/pi;
[U,t,x] = BurgerSolver5(h,k,epsilon,tmax);
%% find closest point to 0
x = 
