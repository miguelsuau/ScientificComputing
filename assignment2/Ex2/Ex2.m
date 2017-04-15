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


%% t
utt = @(x,y) 0*x.*y + 1;
Utt = utt(X,Y);
Uttc = Utt;
Uttc(Iint, Jint) = reshape( poisson9(m)\form_rhs(m, @(x,y) 0*x.*y, utt), m, m );
figure
surf(X,Y,Utt);
figure
surf(X,Y,Uttc);

%% Test case
ut = @(x,y) exp(pi*x).*sin(pi*y)+0.5*(x.*y).^2;
Ut = ut(X,Y);
ft = @(x,y) x.^2+y.^2;

figure(1)
surf(X,Y,Ut);

figure(2)
Utcalc = Ut;
Utcalc(Iint,Jint) = reshape( poisson9(m)\form_rhs(m,ft,ut), m, m );
surf(X,Y,Utcalc);

%% Test case 0
u0 = @(x,y) sin(4*pi*(x+y))+cos(4*pi*x.*y);
U0 = u0(Xint,Yint);

figure(1)
surf(Xint,Yint,U0);

figure(2)
U0calc = reshape( poisson9(m)\form_rhs(m,ft,u0), m, m );
surf(Xint,Yint,U0calc);

%% Test case 1
u1 = @(x,y) x.^2 + y.^2;
U1 = u1(Xint,Yint);

figure(1)
surf(Xint,Yint,U1);

figure(2)
U1calc = reshape( poisson9(m)\form_rhs(m, @(x,y) 0*x.*y, u1), m, m );
surf(Xint,Yint,U1calc);
