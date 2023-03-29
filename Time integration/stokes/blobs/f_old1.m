function Xvp = 2(t,Xv)
% Returns velocity for each blob center

fks = @(t) [0;0];
n = length(Xv)/2;
Xvp = zeros(2*n,1);

% spring forces
% first point
% middle points
% last point
for i = 1:n
    Xvp(2*i-1:2*i) = velocity_regularized(Bp(Xv,i)',Bp(Xv,1),fSpring(Xv,1,n))...
        + velocity_regularized(Bp(Xv,i)',Bp(Xv,n),fSpring(Xv,n,1));
end

% external forcing
for i = 1:n
    Xvp(2*i-1:2*i) = Xvp(2*i-1:2*i) + velocity_regularized(Bp(Xv,i)',Bp(Xv,n),fks(t));
end

end

function [X] = Bp(Xv,i) % returns ith blob position
    X = [Xv(2*i-1);Xv(2*i)];
end

function [F] = fSpring(Xv,i,j) % returns spring force on blob i, from j
    k = 25; % spring constant 
    len_rest = .5; 
    len_curr = norm(Bp(Xv,j)- Bp(Xv,i),2); 
    f_mag = k*(len_curr-len_rest);
    unit_vector = (Bp(Xv,j) - Bp(Xv,i))/len_curr;
    F = f_mag * unit_vector;
end

    