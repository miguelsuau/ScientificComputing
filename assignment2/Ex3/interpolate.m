function [ R ] = interpolate( Rc, m )

mc = (m-1)/2;
Rc = reshape(Rc,mc,mc);
R = zeros(m);

for k=1:mc
    for l=1:mc
        % center
        R(2*k, 2*l) = R(2*k, 2*l) + Rc(k,l);
        % top left
        R(2*k-1, 2*l-1) = R(2*k-1, 2*l-1) + Rc(k,l)/4;
        % top right
        R(2*k-1, 2*l+1) = R(2*k-1, 2*l+1) + Rc(k,l)/4;
        % bottom left
        R(2*k+1, 2*l-1) = R(2*k+1, 2*l-1) + Rc(k,l)/4;
        % bottom right
        R(2*k+1, 2*l+1) = R(2*k+1, 2*l+1) + Rc(k,l)/4;
        % left
        R(2*k, 2*l-1) = R(2*k, 2*l-1) + Rc(k,l)/2;
        % top
        R(2*k-1, 2*l) = R(2*k-1, 2*l) + Rc(k,l)/2;
        % right
        R(2*k, 2*l+1) = R(2*k, 2*l+1) + Rc(k,l)/2;
        % bottom
        R(2*k+1, 2*l) = R(2*k+1, 2*l) + Rc(k,l)/2;
    end
end

end

