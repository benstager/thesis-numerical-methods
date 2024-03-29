function v = velocity_regularized(x,xks,fks)
% Returns velocity vector in R2

v = zeros(2,1);
mu = 1;
n = size(xks,2);
epsilon = .5/10;

for i = 1:n
    v_k = -(fks(:,i))/(4*pi*mu) * ( log(sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)+epsilon)...
            -  (epsilon*(sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)+2*epsilon))/ ...
                ((sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)+epsilon) * sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)) ...
            )...
            + (1/(4*pi*mu))*(dot(fks(:,i),x'-xks(:,i)))*(x'-xks(:,i))*...
            (sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)+2*epsilon)/...
            ((sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)+epsilon)^2 * (sqrt(norm(x'-xks(:,i),2)^2 + epsilon^2)));
    v = v + v_k;

end
