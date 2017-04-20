function [ AU ] = Amult( U, m )
%AMULT for Vcycle

h = 1/(m+1);

% convert to a squared matrix for convenience
U=reshape(U,m+2,m+2);
Iint = 2:m+1;
Jint = 2:m+1;

AU = U;
AU(Iint,Jint) = (1/h^2)*( U(Iint-1,Jint) + U(Iint+1,Jint) + U(Iint,Jint-1) + U(Iint, Jint+1) - 4*U(Iint,Jint) );

% reshape into a column vector -Au
AU = -AU(:);

end

