% exact solution and RHS
u=@(x,y) exp(pi*x).*sin(pi*y)+0.5*(x.*y).^2;
f=@(x,y) x.^2+y.^2;

m=2^6-1;

[X,Y] = meshgrid(linspace(0,1,m+2), linspace(0,1,m+2));

Uall = u(X,Y);
Uall(2:end-1, 2:end-1) = 0;

U = Uall(:);

F = form_rhsV(m,f,u); 

epsilon = 1.0E-10;
omega = 4/5; 

for i=1:100
    R = F(:) + Amult(U,m); % Au = F -> R = F - Au
    fprintf('*** Outer iteration: %3d, rel. resid.: %e\n', ...
        i, norm(R,2)/norm(F,2));
    if(norm(R,2)/norm(F,2) < epsilon)
        break;
    end
    U=Vcycle(U,omega,1000,m,F);
    plotU(m,U);
end

