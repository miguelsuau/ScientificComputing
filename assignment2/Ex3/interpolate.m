function [ R ] = interpolate( Rc, m )

mc = (m-1)/2;
Rc = reshape(Rc, mc+2, mc+2);
R = zeros(m+2,m+2);

% coarsed pts
R(1:2:end, 1:2:end) = Rc;
% inner points
R(2:2:m+1, 2:2:m+1) = 1/4*(Rc(1:end-1, 1:end-1) ...%topleft
    + Rc(2:end,2:end)...%bottomright
    + Rc(1:end-1,2:end)...%top right
    + Rc(2:end,1:end-1)...% bottom left
    );
% add above and below
R(2:2:end-1,1:2:end) = 1/2*(Rc(1:end-1,:) +  Rc(2:end,:));
% left and right
R(1:2:end,2:2:end-1) = 1/2*(Rc(:,1:end-1) +  Rc(:,2:end));

% for k=1:mc
%     for l=1:mc
%         % center
%         R(2*k, 2*l) = R(2*k, 2*l) + Rc(k,l);
%         % top left
%         R(2*k-1, 2*l-1) = R(2*k-1, 2*l-1) + Rc(k,l)/4;
%         % top right
%         R(2*k-1, 2*l+1) = R(2*k-1, 2*l+1) + Rc(k,l)/4;
%         % bottom left
%         R(2*k+1, 2*l-1) = R(2*k+1, 2*l-1) + Rc(k,l)/4;
%         % bottom right
%         R(2*k+1, 2*l+1) = R(2*k+1, 2*l+1) + Rc(k,l)/4;
%         % left
%         R(2*k, 2*l-1) = R(2*k, 2*l-1) + Rc(k,l)/2;
%         % top
%         R(2*k-1, 2*l) = R(2*k-1, 2*l) + Rc(k,l)/2;
%         % right
%         R(2*k, 2*l+1) = R(2*k, 2*l+1) + Rc(k,l)/2;
%         % bottom
%         R(2*k+1, 2*l) = R(2*k+1, 2*l) + Rc(k,l)/2;
%     end
% end

end

