function [p] = pressure_regularized(x,xks,fks)

% Regularized implementation

p = 0;
n = size(xks,2);
epsilon = 0;


for k = 1:n
    pk = (.5*1/pi) * dot(fks(:,k)',x-xks(:,k)')...
        * (norm(x-xks(:,k)',2)^2 + 2*epsilon^2 ...
        + epsilon*sqrt(norm(x-xks(:,k)',2)^2 + epsilon^2)) ...
        / ((sqrt(norm(x-xks(:,k)',2)^2 + epsilon^2))* ...
        (norm(x-xks(:,k)',2)^2 + epsilon^2)^(3/2));
    p = p + pk;
end

