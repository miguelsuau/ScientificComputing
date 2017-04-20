function plotU(m,U)
h=1/(m+1);
x=linspace(0,1,m+2);
y=linspace(0,1,m+2);
[X,Y]=meshgrid(x,y);
surf(X, Y, reshape(U,[m+2,m+2]));
shading interp;
title('Computed solution');
xlabel('x');
ylabel('y');
zlabel('U');
end

