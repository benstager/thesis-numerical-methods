function [u] = velocity_field(n,x)

% n is amount of xks

u = zeros(2,1);
xks = stokes_parameters(n);
fks = forces(xks);
mu = .2;

for i = 1:n
    uk = (-fks(:,i)/(4*pi*mu))*log(norm(x-xks(:,i)',2)) + dot(fks(:,i),x-xks(:,i)')*((x-xks(:,i)')/4*pi*norm(x-xks(:,i)',2))';
    u = u + uk;
end

