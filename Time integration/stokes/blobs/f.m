function Xvp = f(t,Xv)
% Returns velocity for each blob center

fks = @(t) [-.5;0];

n = length(Xv)/2;
Xvp = zeros(2*n,1);
k = 3;

% Finding appropriate way to determine equilibrium position of spring
delta = .07*n;

% left most blob
for i = 1
    Xvp(2*i-1:2*i) = velocity_regularized...
    ([Xv(2*i-1),Xv(2*i)],[Xv(end-1);Xv(end)],fks(t));
end

% middle blobs
for i = 2:n-1
    % contributing forces to left of current cutoff
    f_spring_left = -k*(norm([Xv(2*i-3);Xv(2*i-2)] - [Xv(2*i-1);Xv(2*i)],2)-delta);

    unit_left = [Xv(2*i-3);Xv(2*i-2)] - [Xv(2*i-1);Xv(2*i)] ...
                / norm([Xv(2*i-3);Xv(2*i-2)] - [Xv(2*i-1);Xv(2*i)],2);

    F_spring_left = f_spring_left * unit_left;

    % contributing forces to right of current cutoff
    f_spring_right = -k*(norm([Xv(2*i+1);Xv(2*i+2)] - [Xv(2*i-1);Xv(2*i)],2)-delta);

    unit_right = [Xv(2*i+1);Xv(2*i+2)] - [Xv(2*i-1);Xv(2*i)] ...
                / norm([Xv(2*i+1);Xv(2*i+2)] - [Xv(2*i-1);Xv(2*i)],2);

    F_spring_right = f_spring_right * unit_right;
    
    % total force uses law of superposition
    total_force = F_spring_left + F_spring_right;
    
    % applying force to velocity function
    Xvp(2*i-1:2*i) = velocity_regularized...
    ([Xv(2*i-1),Xv(2*i)],[Xv(1);Xv(2)],total_force);

end

% right most blob
for i = n
    f_spring_left = -k*(norm([Xv(2*i-3);Xv(2*i-2)] - [Xv(2*i-1);Xv(2*i)],2)-delta);
    unit_left = [Xv(2*i-3);Xv(2*i-2)] - [Xv(2*i-1);Xv(2*i)] ...
                / norm([Xv(2*i-3);Xv(2*i-2)] - [Xv(2*i-1);Xv(2*i)],2);
    F_spring_left = f_spring_left * unit_left;

    Xvp(2*i-1:2*i) = velocity_regularized...
    ([Xv(2*i-1),Xv(2*i)],[Xv(1);Xv(2)],F_spring_left);
end


end
