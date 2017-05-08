%% Another nonuniform grid
h = 0.003;
k = h^2;
tmax = 1.6037/pi;
xspan = [-1 1];
epsilon = 0.01/pi;

nonunifpoints = 35;

[Unu,tnu,xnu] = BurgerSolverNonUniform(nonunifpoints,k,epsilon,tmax);

idx = floor(length(xnu)/2); % mid point x = 0
dxnu = (Unu(idx+1,end)-Unu(idx-1,end))/(xnu(idx+1)-xnu(idx-1));
fprintf('Non uniform grid -> dx_nu = %f\n', dxnu);

figure(1);
% Show the solution
subplot(2,1,1);
plot(xnu,Unu(:,end));
title('Calculated solution on the non-uniform grid');
axis([-1 1 -1.1 1.1]); grid on;
xlabel('x'); ylabel('u');
% Show the grid
subplot(2,1,2);
plot(xnu,zeros(length(xnu)),'o-');
title(sprintf('Visualization of the non-uniform grid for %d points', nonunifpoints));

%% True solution
h = 0.003;
k = h^2;
epsilon = 0.01;
tmax = 1;

nonunifpoints = 35;

[Unu,tnu,xnu] = BurgerSolverNonUniform(nonunifpoints,k,epsilon,tmax);

x = nonuniformgrid(nonunifpoints);
t = linspace(0,tmax,tmax/k);
[X,T] = ndgrid(x,t);
Utrue = boundaryFun(X,T,epsilon);

figure(1);
plot(x,Unu(:,end),'ro-',x,Utrue(:,1),'b-');
