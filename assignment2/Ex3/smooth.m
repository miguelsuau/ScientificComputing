function Unew=smooth(U,omega,m,F)
% return m+2 by m+2
    h = 1/(m+1);
    U = reshape(U,m+2,m+2);
    F = reshape(F,m+2,m+2);
    Iint = 2:m+1;
    Jint = 2:m+1;
    U(Iint,Jint) = (1-omega)*U(Iint,Jint) + omega/4*( U(Iint-1,Jint) + U(Iint+1,Jint) + U(Iint,Jint-1) + U(Iint, Jint+1) - h^2*F(Iint,Jint));
    Unew = U(:);
end