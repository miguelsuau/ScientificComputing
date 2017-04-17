%% Setup
close all;

a = 0; 
b = 1; 
m = 20;
h = (b-a)/(m+1);
x = linspace(a,b,m+2);   % grid points x including boundaries
y = linspace(a,b,m+2);   % grid points y including boundaries


[X,Y] = meshgrid(x,y);      % 2d arrays of x,y values

Iint = 2:m+1;              % indices of interior points in x
Jint = 2:m+1;              % indices of interior points in y
Xint = X(Iint,Jint);       % interior points
Yint = Y(Iint,Jint);

fdummy = @(x,y) 0*x.*y;

%% Dummy test
utt = @(x,y) 0*x.*y + 1;
Utt = utt(X,Y);
Uttc = Utt;
Uttc(Iint, Jint) = reshape( poisson9(m)\form_rhs(m, @(x,y) 0*x.*y, utt), m, m );
figure
surf(X,Y,Utt);
figure
surf(X,Y,Uttc);
figure
semilogy(max( abs( Utt-Uttc ) ))

%% Test case 0
u0 = @(x,y) sin(4*pi*(x+y))+cos(4*pi*x.*y);
U0 = u0(X,Y);

figure(1)
surf(X,Y,U0);

figure(2)
U0calc = U0;
U0calc(Iint,Jint) = reshape( poisson9(m)\form_rhs(m, @(x,y) 0*x.*y, u0), m, m );
surf(X,Y,U0calc);
%TODO global error
figure(3)
semilogy(max( abs( U0-U0calc ) ))

%% Test case 1 (unfinished)
u1 = @(x,y) x.^2 + y.^2;
U1 = u1(Xint,Yint);

figure(1)
surf(Xint,Yint,U1);

figure(2)
U1calc = reshape( poisson9(m)\form_rhs(m, @(x,y) x.^2+y.^2, u1), m, m );
surf(Xint,Yint,U1calc);

%% Test case 2 (unfinished)
u2 = @(x,y) sin(2*pi*abs(x - y).^(2.5));
U2 = u2(X,Y);

figure(1)
surf(X,Y,U2);

figure(2)
U2calc = U2;
U2calc(Iint,Jint) = reshape( poisson9(m)\form_rhs(m, u2, u2), m, m );
surf(X,Y,U2calc);
