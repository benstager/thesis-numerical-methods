function [p] = pressure(x,xks,fks)

% suppose n is amount of xks, suppose function only returns pressure at one
% sample point

p = 0;
n = size(xks,2);

for k = 1:n
    pk = dot(fks(:,k)',(x-xks(:,k)'))/(norm(x-xks(:,k)',2)^2);
    p = p + pk;
end

