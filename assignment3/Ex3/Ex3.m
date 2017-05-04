%% Non-linear advection-diffusion equation
h = 0.1;
k = h^2/6;
epsilon = 0.5;
tmax = 1;
U = BurgerSolver(@boundaryFun,h,k,epsilon,tmax);
mesh(U);
figure
%% true solution
x = linspace(-1,1,100);
t = linspace(0,1,100);
[X,T] = meshgrid(x,t);
U = boundaryFun(X,T,epsilon);
mesh(X,T,U)