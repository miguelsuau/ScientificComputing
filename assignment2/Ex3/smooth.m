function Unew=smooth(U,omega,m,F)
    h = 1/(m+1);
    UJ = U;
    UJ(2:end-1,2:end-1) = 0.25*(U(1:end-2,2:end-1)+U(3:end,2:end-1)+ ...
              U(2:end-1,1:end-2)+U(2:end-1,3:end)-h^2.*F(2:end-1,2:end-1));
    Unew = U + omega.*(UJ-U);
end