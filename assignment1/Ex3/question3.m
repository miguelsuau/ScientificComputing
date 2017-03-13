% Calculate our method's Butcher's tableau
format rat
f = @(x) [  x(1)+x(2)+x(3)-1;
            1/4*x(2)+x(3)-1/2;
            1/16*x(2)+x(3)-1/3;
            1/4*x(3)*x(6)-1/6;
            x(4)-1/4;
            x(5)+x(6)-1;
         ];
x = fsolve(f,zeros(6,1));
b1 = x(1); b2 = x(2); b3 = x(3); a21 = x(4); a31 = x(5); a32 = x(6);

butcher.AT = [0 0 0; a21 0 0; a31 a32 0]';
butcher.b  = [b1 b2 b3]';
butcher.c  = [0 1/3 1]';
butcher.d  = [b1 b2 b3]' - [1/8 1/2 3/8]'; % b - bhat
butcher.stages = 3;

% Test equation
addpath('../Ex1'); % if ran from Ex3 directory 
tspan = [0; 10];
n = 20; % change as needed
x0 = 1;
lambda = -1;

[T1,X1,Err1] = ExplicitRungeKutta(@TestEquation,tspan,x0,n,butcher,lambda);

% Analytical solution
T = linspace(tspan(1),tspan(2),n)';
X = exp(lambda*T);

% Plot method vs 'true' solution
plot(T1,X1(:),'-g','LineWidth',2.5);
hold on
plot(T,X(:),'--r','LineWidth',2.5);
legend('Analytical solution', 'Designed Runge-Kutta')
hold off

% TODO: Discuss how the local error (as a function of step size) should be
% plotted - i.e. how many graphs? subplots?
alpha = -5:0.01:5;
beta = -5:0.01:5;
A = butcher.AT';
b = butcher.b;
c = butcher.c;
d = butcher.d;
nreal = length(alpha);
nimag = length(beta);
I = eye(size(A));
e = ones(size(A,1),1);

for kreal = 1:nreal
    for kimag = 1:nimag
    z = alpha(kreal) + 1i*beta(kimag);
    tmp = (I-z*A)\e;
    R = 1 + z*b'*tmp;
    Ehat = z*d'*tmp;
    f = exp(z);
    E = R-f;
    EhatmE = Ehat-E;
    absR(kimag,kreal) = abs(R);
    absEhatmE(kimag,kreal) = abs(EhatmE);
    absEhat(kimag,kreal) = abs(Ehat);
    absE(kimag,kreal) = abs(E);
    absF(kimag,kreal) = abs(f);
    end
end

figure
fs = 14;
subplot(221)
imagesc(alpha,beta,absR,[0 1]);
grid on
colorbar
axis image
axis xy
xlabel('real','fontsize',fs);
ylabel('imag','fontsize',fs);
title('|R(z)|','fontsize',fs)
