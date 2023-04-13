function v = velocity_regularized_3d(x,xks,fks)
% Returns velocity vector in R2

v = zeros(3,1);
mu = 1;
n = size(xks,2);
epsilon = .5/10;
r = @(x,y) norm(x-y,2);

for i = 1:n
    v_k = (1/mu)*(r(x',xks(:,i))^2+2*epsilon^2)...
    /(8*pi*(r(x',xks(:,i))^2+epsilon^2)^(3/2))*fks(:,i)...
    + (dot(fks(:,i),(x'-xks(:,i))))*(x'-xks(:,i))*...
    1/(8*pi*(r(x',xks(:,i))^2+epsilon^2))^(3/2);
    
    v = v + v_k;
end
