function U = NewtonsMethod(FunJac,Uinit,h,tol,maxit)
U = Uinit;
[G,J] = FunJac(U,h);
k = 0;
while ((k < maxit) && (norm(G,'inf') > tol))
    k = k+1;
    U(2:end-1) = U(2:end-1) - J\G;
    [G,J] = FunJac(U,h);
end