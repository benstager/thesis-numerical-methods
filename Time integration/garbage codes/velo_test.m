function v = velo_test(x,xks,fks)

v = zeros(2,1);
mu = 1;
n = size(xks,2);


for i = 1:n
    v_k = (-fks(:,i)/(4*pi*mu))*log(norm(x'-xks(:,i),2))...
    +(dot(fks(:,i),x'-xks(:,i)))*(x'-xks(:,i))/(4*pi*mu*norm(x'-xks(:,i),2)^2);
    v = v + v_k;
end
function [u] = velocity(x,xks,fks)

% Function returns velocity for certain 'x' coordinate in R2

u = zeros(2,1);
mu = 1;
n = size(xks,2);

for i = 1:n
    uk = (-fks(:,i)/(4*pi*mu))*log(norm(x-xks(:,i)',2)) + dot(fks(:,i),x-xks(:,i)')*((x-xks(:,i)')/4*pi*norm(x-xks(:,i)',2))';
    u = u + uk;
end
