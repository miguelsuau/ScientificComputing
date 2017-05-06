%% Non-linear advection-diffusion equation
h = 0.1;
k = h^2/6;
epsilon = 0.5;
tmax = 1;
U1 = BurgerSolver(@boundaryFun,h,k,epsilon,tmax);
mesh(U1);
figure
%% true solution
x = linspace(-1,1,100);
t = linspace(0,1,100);
[X,T] = ndgrid(x,t);
U2 = boundaryFun(X,T,epsilon);
mesh(X,T,U2)
%% 3.3
h = 0.1;
k = h^2;
tmax = 1.6037/pi;
epsilon = 0.01/pi;
U3 = BurgerSolver2(h,k,epsilon,tmax);