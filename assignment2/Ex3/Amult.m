function [ AU ] = Amult( U, m )
%AMULT Summary of this function goes here
%   Detailed explanation goes here

h = 1/(m+1);

AU = zeros(m,m);

for k=1:m
    for l=1:m
        % calculate the (AU)_{ij} based on 5-point L
        % convert indexing from (i,j) to (i-1)*m+j since U is a vector
        % or i*m-m+j or i*m-(m-j)
        Us = -4*U(k*m-m+l); % (i,j) always present
        % check if the index of i is either 0
        if (k ~= 1)
            Us = Us + U((k-1)*m-m+l); % include the top point
        end
        % or mk
        if (k ~= m)
            Us = Us + U((k+1)*m-m+l); % include the bottom point
        end
        % check if the index of j is either 0
        if (l ~= 1)
            Us = Us + U(k*m-m+(l-1));
        end
        % or ml
        if (l ~= m)
            Us = Us + U(k*m-m+(l+1));
        end
        % finalize
        AU(k,l) = (1/h^2)*Us;
    end
end

end

