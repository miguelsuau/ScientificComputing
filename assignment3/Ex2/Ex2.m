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

%% Visualize `u`
[m,n] = size(U);
maskt = 1:501;
maskx = 1:n;

figure(1)
subplot(2,2,1)
contour(X(maskt,maskx),T(maskt,maskx),U(maskt,maskx)), colorbar
xlabel('x'), ylabel('t')

subplot(2,2,2)
mesh(X(maskt,maskx),T(maskt,maskx),U(maskt,maskx))
xlabel('x'), ylabel('t'), zlabel('u')

subplot(2,2,3)
plot(T(maskt,1),U(maskt,1:1:5))
xlabel('t'), ylabel('U(t,x=1,2,3,4,5)')

subplot(2,2,4)
plot(X(1,:), U(1:1:5, :))
xlabel('x'), ylabel('U(x,t=1,2,3,4,5)')

%% FTBS
Uftbs = ftbs(u, cr, 40, 100);

[m,n] = size(Uftbs);
maskt = 1:501;
maskx = 1:n;

figure(2)
subplot(2,2,1)
contour(X(maskt,maskx),T(maskt,maskx),Uftbs(maskt,maskx)), colorbar
xlabel('x'), ylabel('t')

subplot(2,2,2)
mesh(X(maskt,maskx),T(maskt,maskx),Uftbs(maskt,maskx))
xlabel('x'), ylabel('t'), zlabel('u')

subplot(2,2,3)
plot(T(maskt,1),Uftbs(maskt,1:1:5))
xlabel('t'), ylabel('U(t,x=1,2,3,4,5)')

subplot(2,2,4)
plot(linspace(-1,1,n), Uftbs(1:1:5, :))
xlabel('x'), ylabel('U(x,t=1,2,3,4,5)')

% Show comparison
figure(1); figure(2);

figure(4)
xx = linspace(-1,1,n);
Ut = U(:,1:end-1);
subplot(3,1,1)
plot(xx, Uftbs(1, :), 'bo-', xx, Ut(1,:), 'rx-');
subplot(3,1,2)
plot(xx, Uftbs(1001, :), 'bo-', xx, Ut(1001,:), 'rx-');
subplot(3,1,3)
plot(xx, Uftbs(2001, :), 'bo-', xx, Ut(2001,:), 'rx-');

%% Upwind
ShowMovingSolution = false;
nWavePeriods = 40;

a = 1/2;        % speed
f = 1/2;        % frequency
T = f^(-1);     % period
lambda = 1;     % wave length
dx = 1/100;
dt = 16/1000;
Cr = a*dx/dt;
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

%% Analysis of dispersion and diffusion
figure(7);
g40 = g(0:40, 0.01);
plot(0:40, abs(g40), 'o-');
title('Upwind method -- diffusion with $\xi = 1,2,\dots,40$; $\Delta x = \frac{1}{100}$; $\nu = \frac{8}{10}$', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$\Re \big( g(\xi) \big)$', 'Interpreter', 'latex', 'FontSize', 16)

figure(8);
plot(0:40, angle(g40), 'x-');
title('Upwind method -- dispersion with $\xi = 1,2,\dots,40$; $\Delta x = \frac{1}{100}$; $\nu = \frac{8}{10}$', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$\Im \big( g(\xi) \big)$', 'Interpreter', 'latex', 'FontSize', 16)

%% Analysis of the errors
nWavePeriods = 40;

a = 1/2;        % speed
f = 1/2;        % frequency
T = f^(-1);     % period
lambda = 1;     % wave length
N = 200;        % number of points - 1
dx = 1/100;
dt = 16/1000;
xmin = -1;
xmax = 1;
t = 0;
tmax = nWavePeriods*T; % tmax = 80 s
steps = tmax/dt;

x = xmin-dx:dx:xmax+dx;


[Uup, Utrue] = upwind(u, tmax);

