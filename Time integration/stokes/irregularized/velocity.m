function [u] = velocity(x,xks,fks)

% Function returns velocity for certain 'x' coordinate in R2

u = zeros(2,1);
mu = .2;
n = size(xks,2);

for i = 1:n
    uk = (-fks(:,i)/(4*pi*mu))*log(norm(x-xks(:,i)',2)) + dot(fks(:,i),x-xks(:,i)')*((x-xks(:,i)')/4*pi*norm(x-xks(:,i)',2))';
    u = u + uk;
end

