function Unew=Vcycle(U,omega,nsmooth,m,F)
% Approximately solve: A*U = F
h=1.0/(m+1);
l2m=log2(m+1);
assert(l2m==round(l2m));
assert(length(U)==m*m);

if(m==1)
    % if we are at the coarsest level
    % TODO: solve the only remaining equation directly!
    fprintf('Returning U: %f and F: %f \n', U, F);
    Unew = Amult(U,m)\F;
else
    fprintf('Running recursive part m: %d \n', m);
    % 1. TODO: pre-smooth the error
    %    perform <nsmooth> Jacobi iterations
    Unew = U;
    for numSmooth=1:nsmooth
        Unew = smooth(Unew,omega,m,F);
    end
    % 2. TODO: calculate the residual
    R = F(:) + Amult(Unew,m); % Amult returns -Au, so +
    % 3. TODO: coarsen the residual
    Rcoarse = coarsen(R,m);
    % 4. recurse to Vcycle on a coarser grid
    mc=(m-1)/2;
    Ecoarse=Vcycle(zeros(mc*mc,1),omega,nsmooth,mc,-Rcoarse);
    % 5. TODO: interpolate the error
    E = interpolate(Ecoarse,m);
    % 6. TODO: update the solution given the interpolated error
    Unew = Unew + E(:);
    % 7. TODO: post-smooth the error
    %    perform <nsmooth> Jacobi iterations
    for numSmooth=1:nsmooth
        Unew = smooth(Unew,omega,m,F);
    end
end

end