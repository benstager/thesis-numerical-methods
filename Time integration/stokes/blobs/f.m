function Xvp = f(t,Xv)
% Returns velocity for each blob center

fks = @(t) [1;0];
n = length(Xv)/2;
Xvp = zeros(2*n,1);

for i = 1:n
    Xvp(2*i-1:2*i) = velocity_regularized([Xv(2*i-1),Xv(2*i)],[0;0],fks(t));
end

end
