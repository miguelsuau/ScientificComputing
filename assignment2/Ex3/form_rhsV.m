function [ F ] = form_rhsV( m, f, u )
%FORM_RHSV Summary of this function goes here
%   Detailed explanation goes here

h = 1/(m+1);
[X,Y] = meshgrid(linspace(0,1,m+2), linspace(0,1,m+2));
Iint = 2:m+1;
Jint = 2:m+1;
Xint = X(Iint,Jint);
Yint = Y(Iint,Jint);

F = u(X,Y);
F(Iint,Jint) = f(Xint,Yint);

F=F(:);

end

