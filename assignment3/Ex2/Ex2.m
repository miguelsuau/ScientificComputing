%% von Neumann analysis of the Upwind method
v = 0.8; % Courant number
g = @(xi,h) 1 - v + v*exp(-1i*xi*h);
figure(1);
clf;
rectangle('Position',[-1 -1 2 2],'Curvature',[1 1], 'EdgeColor', 'r', 'LineStyle', '-.', 'LineWidth', 1.2);
hold on;
plot(g(-1000:0.01:1000,0.01), 'LineWidth', 2); % verify that values are bound within a unit circle
hold off;
axis equal;
grid on;
title('Upwind method -- von Neumann analysis with $\nu = 0.8$', 'Interpreter', 'latex', 'FontSize', 16);
legend({'g(\xi)'}, 'FontSize', 14);

%% Define the problem
f = @(x) sin(2*pi*x);
cr = 0.8;
u = @(x,t) f(x-0.5*t);

a=1/2;
wave_periods = 40;
n_wave_length = 100;
x = linspace(-1,1,100+2); % \Delta x = 1/100
t = linspace(0,a^(-1)*wave_periods,a^(-1)*n_wave_length*wave_periods);

[X,T] = meshgrid(x,t);
U = u(X,T);
U(:,1) = u(1,t); % at fixed x
U(:,end) = u(-1,t); % at fixed x
for k=1:(2*n_wave_length):length(t)
    U(k,:) = u(x,0); % periodic BCs
end

%% Animation of the solution at a const speed
figure
h = animatedline;
axis([-1 1 -1 1])
legend('');

for k = 1:length(t)
    h.DisplayName = num2str(round(T(k,1), 2));
    clearpoints(h)
    addpoints(h,X(k,:),U(k,:))
    drawnow
end

%% Upwind
ShowMovingSolution = false;
nWavePeriods = 40;

a = 1/2;        % speed
f = 1/2;        % frequency
T = f^(-1);     % period
lambda = 1;     % wave length
dx = 1/100;
dt = 16/1000;
Cr = a*(dt/dx);
xmin = -1;
xmax = 1;
t = 0;
tmax = nWavePeriods*T; % tmax = 80 s
steps = tmax/dt;

x = xmin-dx:dx:xmax+dx;


[Uup, Utrue] = upwind( u, a, dx, dt, xmin, xmax, tmax );
[nsteps, nx] = size(Uup);

% Indices of 1-40 wave periods
uIdxRange = 1:(steps/nWavePeriods):steps+1;
[um,un] = size(Uup);
xx = linspace(xmin,xmax,un);

figure(5);
plot(xx, Uup(uIdxRange(end), :), 'ro',...
     xx, Utrue(1,   :), 'b-',...
    'LineWidth', 1.2)
legend({'$u(x,80)$ upwind', '$u(x,80)$'}, 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast')
xlabel('x', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('u', 'Interpreter', 'latex', 'FontSize', 16);
grid on;


if ShowMovingSolution
    figure(6); clf;
    % Show step to t=3.6s
    for n=1:(3.6/dt)
        t = t+dt;
        figure(7);
        plot(x(2:end-1),Uup(n,2:end-1),'bx-'); hold on
        plot(x(2:end-1),Utrue(n,2:end-1),'g-','LineWidth',2); hold off
        axis([-1 1 -1 1]);
        title(sprintf('Time = %2.2f s', t))
        pause(0.0001);
    end
end

%% Analysis of the errors
nWavePeriods = 40;

a = 1/2;        % speed
f = 1/2;        % frequency
T = f^(-1);     % period
lambda = 1;     % wave length
dx = 1/100;
dt = 16/1000;
Cr = 0.8; % fixed
xmin = -1;
xmax = 1;
t = 0;
tmax = nWavePeriods*T; % tmax = 80 s
steps = tmax/dt;
err = [];
dxs = []; dts = [];
g = @(phi) 1 - Cr + Cr*exp(1i*phi); % phi = xi*dx
xi = 2*pi/1; % 2*pi/L, L = 1

fprintf('Running the convergence check loop...\n');
for vdt = [1/100 1/150 1/200 1/250 1/300]
    vdx = a*vdt/Cr; % Cr = a*(dt/dx) => dx = a*dt/Cr
    fprintf('Using Cr = %f, dt = %f, dx = %f\n', a*vdt/vdx, vdt, vdx);
    [Uup, Utrue] = upwind( u, a, vdx, vdt, xmin, xmax, tmax );
    err = [err max(abs(Uup(end,:)-Utrue(end,:)))];
    dxs = [dxs vdx]; dts = [dts vdt];
end
figure(1); clf;
subplot(1,2,1);
loglog(dts,dts,dts,err,'o-');
legend({'$\mathcal{O}(\Delta t)$', 'Upwind'},'Interpreter','latex','location','southeast','FontSize',16)
set(gca, 'FontSize', 14);
subplot(1,2,2);
loglog(dxs,dxs,dxs,err,'o-');
legend({'$\mathcal{O}(\Delta x)$', 'Upwind'},'Interpreter','latex','location','southeast','FontSize',16)
set(gca, 'FontSize', 14);

%% Finding the phase difference
[c,lag] = xcorr(Uup(end,:), Utrue(1,:));
[maxC,cIdx]=max(c);
fprintf('Phase difference is %e\n', lag(cIdx));

