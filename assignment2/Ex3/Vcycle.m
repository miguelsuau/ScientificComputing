function Unew=Vcycle(Ufull,omega,nsmooth,m,F)
% Approximately solve: A*U = F
h=1.0/(m+1);
l2m=log2(m+1);
assert(l2m==round(l2m));

assert(length(Ufull)==(m+2)^2);

if(m==1)
    U = reshape(Ufull, 3, 3);
    Uval = 1/4 * (U(1,2) + U(2,1) + U(2,3) + U(3,2) - h^2*F(5));
    Unew = U; 
    Unew(2,2) = Uval;
    Unew = Unew(:);
else
    %fprintf('Running the recursive part m: %d \n', m);

    Unew = Ufull;
    for numSmooth=1:nsmooth
        Unew = smooth(Unew,omega,m,F(:));
    end
    
    R = F(:) + Amult(Unew(:),m); % Amult returns -Au, so +
    
    Rcoarse = coarsen(R,m);
    
    mc=(m-1)/2;
    Ecoarse=Vcycle(zeros((mc+2)^2,1),omega,nsmooth,mc,-Rcoarse);
    
    %fprintf('Returning from the recursive part m: %d \n', m);
    
    E = interpolate(Ecoarse(:),m);
    E = reshape(E, m+2,m+2);
    Err = zeros(size(E));
    Err(2:end-1,2:end-1) = E(2:end-1,2:end-1);
    
    Unew = Unew(:) + Err(:);

    for numSmooth=1:nsmooth
        Unew = smooth(Unew,omega,m,F);
    end
    
    Unew = Unew(:);
end

end