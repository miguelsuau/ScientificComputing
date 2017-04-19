function plotU(m,U)
h=1/(m+1);
x=linspace(1/h,1-1/h,m);
y=linspace(1/h,1-1/h,m);
[X,Y]=meshgrid(x,y);
surf(X, Y, reshape(U,[m,m])');
shading interp;
title('Computed solution');
xlabel('x');
ylabel('y');
zlabel('U');
end

