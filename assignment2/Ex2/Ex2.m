close all;
%% Setup
a = 0; 
b = 1; 
m = 30;
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
utrue = @(x,y) exp(x+y/2);
f = @(x,y) 1.25 * exp(x+y/2);
Utrue = utrue(X,Y);

RHS = form_rhs(m,f,utrue);
A = poisson9(m);
U = Utrue;
U(Iint,Jint) = reshape(A\RHS(:), m, m);

figure(1);
surf(X,Y,Utrue);

figure(2);
surf(X, Y, U);

figure(3)
semilogy(max( abs( Utrue-U ) ))

%% Test case 0
u0 = @(x,y) sin(4*pi*(x+y))+cos(4*pi*x.*y);
U0 = u0(X,Y);
% given by matlab simplify(laplacian() or diff(...,2)+diff(...,2))
f0 = @(X,Y) -16 * (pi ^ 2) * (...
    cos((4 * pi * X .* Y)) .* (X .^ 2)...
    + cos((4 * pi * X .* Y)) .* (Y .^ 2)...
    + 2 * sin((4 * pi * (X + Y)))...
);

RHS0 = form_rhs(m, f0, u0);
A = poisson9(m);
U0calc = U0;
U0calc(Iint,Jint) = reshape( A\RHS0(:), m, m );

figure(1)
surf(X,Y,U0);

figure(2)
surf(X,Y,U0calc);

%% Test case 0 - global error

Gerr = [];
Grange = [30 100 500 1000];
for k=Grange
    m = k;
    h = 1/(m+1);
    x = linspace(0,1,m+2);
    y = linspace(0,1,m+2);
    [X,Y] = meshgrid(x,y);
    Iint = 2:m+1; 
    Jint = 2:m+1;
    Xint = X(Iint,Jint);
    Yint = Y(Iint,Jint);
    
    U0 = u0(X,Y);
    RHS0 = form_rhs(m, f0, u0);
    A = poisson9(m);
    U0calc = U0;
    U0calc(Iint,Jint) = reshape( A\RHS0(:), m, m );
    
    Gerr = [Gerr(:); max(max( abs( U0(Iint,Jint)-U0calc(Iint,Jint) ) ))];
    fprintf('Finished calculating for m=%d\n', m)
end
figure(3)
hs = 1./(10*(Grange+ones(size(Grange)))); % h = 1/(m+1)
loglog(hs, Gerr, hs, hs.^(4+1))
fprintf('Test case 0 done!\n')

%% Test case 1 (unfinished)
u1 = @(x,y) x.^2 + y.^2;
U1 = u1(Xint,Yint);

RHS1 = form_rhs(m, @(x,y) 4*x.^0.*y.^0, u1); % f(x,y) = 4
U1calc = reshape( poisson9(m)\RHS1(:), m, m );

figure(1)
surf(Xint,Yint,U1);

figure(2)
surf(Xint,Yint,U1calc);

figure(3)
semilogy(max( abs( U1-U1calc ) ))

%% Test case 2 (unfinished)
u2 = @(x,y) sin(2*pi*abs(x - y).^(2.5));
U2 = u2(X,Y);
f2 = @(x1,y1) 5*pi*abs(x1 - y1).^(1/2).*sign(x1 - y1).^2.*(...
    3*cos(2*pi*abs(x1 - y1).^(5/2))...
    - 10*pi*abs(x1 - y1).^(5/2).*sin(2*pi*abs(x1 - y1).^(5/2))...
);
U2calc = U2;
U2calc(Iint,Jint) = reshape( poisson9(m)\form_rhs(m, f2, u2), m, m );

figure(1)
surf(X,Y,U2);

figure(2)
surf(X,Y,U2calc);

figure(3)
semilogy(max( abs( U2(Iint,Jint)-U2calc(Iint,Jint) ) ))
