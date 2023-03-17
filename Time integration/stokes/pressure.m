function [p] = pressure(x,xks,fks)

% p function takes in a bivariate 'x' point in R2, k amount of forces, k amount
% of Bivariate points where forces occur, and what those forces are

% Returns a scalar value p, can be looped to find p over discretized square domain

p = 0;
n = size(xks,2);

for k = 1:n
    pk = dot(fks(:,k)',(x-xks(:,k)'))/(norm(x-xks(:,k)',2)^2);
    p = p + pk;
end

