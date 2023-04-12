function Xvp = f(t,Xv)
% Returns velocity for each blob center

fks_right = @(t) [0;10*sin(2*t)];
fks_left = @(t) [0;0];
n = length(Xv)/2;
Xvp = zeros(2*n,1);

% spring forces
% first point, external forcing and last point
for i = 1:n
    Xvp(2*i-1:2*i) = Xvp(2*i-1:2*i) + ...
    velocity_regularized(Bp(Xv,i)',Bp(Xv,n),fks_right(t)+fSpring(Xv,n,n-1))+...
        velocity_regularized(Bp(Xv,i)',Bp(Xv,1),fks_left(t)+ fSpring(Xv,1,2));
        
end

% middle points
for i = 2:n-1
    % contribution of forces from ith blob
    for j = 1:n
        % at j locations
        Xvp(2*j-1:2*j) = Xvp(2*j-1:2*j) + ...
        velocity_regularized(Bp(Xv,j)',Bp(Xv,i), ...
            fSpring(Xv,i,i-1)+fSpring(Xv,i,i+1));
    end
end

end

    