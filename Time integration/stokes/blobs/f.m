function Xvp = f(t,Xv)
% Returns velocity for each blob center

fks = @(t) [.5*(sin(t) +1);.25];

n = length(Xv)/2;
Xvp = zeros(2*n,1);
k = 10;
delta = [.17*(n/2);0];

for i = 1
    Xvp(2*i-1:2*i) = velocity_regularized...
    ([Xv(2*i-1),Xv(2*i)],[Xv(end-1);Xv(end)],fks(t));
end

for i = 2:n-1
    f_spring_left = -k*(norm([Xv(2*i+1);Xv(2*i+2)] - [Xv(2*i-1);Xv(2*i)],2));
    unit_left = [Xv(2*i+1);Xv(2*i+2)] - [Xv(2*i-1);Xv(2*i)] ...
                / norm([Xv(2*i+1);Xv(2*i+2)] - [Xv(2*i-1);Xv(2*i)],2);
    F_spring_left = f_spring_left * unit_left;
    %norm([Xv(2*i-1);Xv(2*i)] - [Xv(2*i-3);Xv(2*-2)],2);
    f_spring_right = -k*(norm([Xv(2*i+1);Xv(2*i+2)] - [Xv(2*i-1);Xv(2*i)],2));
    unit_right = [Xv(2*i+1);Xv(2*i+2)] - [Xv(2*i-1);Xv(2*i)] ...
                / norm([Xv(2*i+1);Xv(2*i+2)] - [Xv(2*i-1);Xv(2*i)],2);
    F_spring_right = f_spring_right * unit_right;

end

for i = n
end


end
