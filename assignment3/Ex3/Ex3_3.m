u = @(x,t,epsi) -tanh((x+0.5-t)/(2*epsi)) + 1;
x = linspace(-1,1,500);
t = linspace(0,2*(1.6037/pi),1000); % 1.6037/pi ~= 0.5105
[X,T] = meshgrid(x,t);
U = u(X,T,0.01/pi);
figure(1);
clf;
mesh(X,T,U); hold on;
dummypts = 0:0.01:2;
scatter3(0*ones(size(dummypts)),1.6037/pi*ones(size(dummypts)), dummypts); hold off;
xlabel('x'); ylabel('t'); zlabel('u');

