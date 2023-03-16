function [f] = velocity_functions(xks,fks)
% THIS IS WRONG- just messing around, how can we solve this dynamical
% system

u = zeros(2,1);
mu = .2;
n = size(xks,2);

for i = 1
   f = @(x,y) (-fks(1,i)/(4*pi*mu))*log(norm([x,y]-xks(:,i)',2)) + dot(fks(1,i),[x,y]-xks(:,i)')*(([x,y]-xks(:,i)')/4*pi*norm([x,y]-xks(:,i)',2))';
end

for i = 1
   g = @(x,y) (-fks(2,i)/(4*pi*mu))*log(norm([x,y]-xks(:,i)',2)) + dot(fks(2,i),[x,y]-xks(:,i)')*(([x,y]-xks(:,i)')/4*pi*norm([x,y]-xks(:,i)',2))';
end



