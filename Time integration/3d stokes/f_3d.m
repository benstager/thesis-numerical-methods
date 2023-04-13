function Xvp = f_3d(t,Xv)
% Returns velocity for each blob center

fks_right = @(t) [0;0;10*sin(2*t)];
fks_left = @(t) [0;0;0];
n = length(Xv)/3;
Xvp = zeros(3*n,1);

% spring forces
% first point, external forcing and last point
for i = 1:n
    Xvp(3*i-2:3*i) = Xvp(3*i-2:3*i) + ...
    velocity_regularized_3d(Bp_3d(Xv,i)',Bp_3d(Xv,n),fks_right(t)+fSpring_3d(Xv,n,n-1))+...
        velocity_regularized_3d(Bp_3d(Xv,i)',Bp_3d(Xv,1),fks_left(t)+ fSpring_3d(Xv,1,2));
end

% middle points
for i = 2:n-1
    % contribution of forces from ith blob
    for j = 1:n
        % at j locations
        Xvp(3*i-2:3*i) = Xvp(3*i-2:3*i) + ...
        velocity_regularized_3d(Bp_3d(Xv,j)',Bp_3d(Xv,i), ...
            fSpring_3d(Xv,i,i-1)+fSpring_3d(Xv,i,i+1));
    end
end

end

    