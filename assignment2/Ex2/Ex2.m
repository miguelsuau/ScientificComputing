close all;

a = 0; 
b = 1; 
m = 20;
h = (b-a)/(m+1);
x = linspace(a,b,m+2);   % grid points x including boundaries
y = linspace(a,b,m+2);   % grid points y including boundaries


[X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
X = X';                     % transpose so that X(i,j),Y(i,j) are
Y = Y';                     % coordinates of (i,j) point

Iint = 2:m+1;              % indices of interior points in x
Jint = 2:m+1;              % indices of interior points in y
Xint = X(Iint,Jint);       % interior points
Yint = Y(Iint,Jint);


u0 = @(x,y) exp(pi*x).*sin(pi*y)+0.5*(x.*y).^2;
U0 = u0(Xint,Yint);
f0 = @(x,y) x.^2+y.^2;

figure(1)
surf(Xint,Yint,U0);

figure(2)
U0calc = reshape( poisson9(m)\form_rhs(m,f0,u0), m, m );
surf(Xint,Yint,U0calc');