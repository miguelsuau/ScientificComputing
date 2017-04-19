function [ AU ] = Amult( U, m )
%AMULT Summary of this function goes here
%   Detailed explanation goes here

h = 1/(m+1);

% convert to a squared matrix for convenience
U=reshape(U,m,m);

AU=4*U;
AU(1:m-1,:) = AU(1:m-1,:) - U(2:m,:);
AU(2:m,:) = AU(2:m,:) - U(1:m-1,:);
AU(:,2:m) = AU(:,2:m) - U(:,1:m-1);
AU(:,1:m-1) = AU(:,1:m-1) - U(:,2:m);

% reshape into a column vector -Au
AU = AU(:)/h^2;

end

