function v = velo_regularized_test(x,xks,fks)
%VELO Summary of this function goes here
%   Detailed explanation goes here

v = zeros(2,1);
mu = 1;
n = size(xks,2);
epsilon = .5;

for i = 1:n
    v_k = -(fks(:,i))/(4*pi*mu) * ( log(sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)+epsilon)...
            -  (epsilon*(sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)+2*epsilon))/ ...
                ((sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)+epsilon) * sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)) ...
            )...
            + (1/(4*pi*mu))*(dot(fks(:,i),x'-xks(:,i)))*(x'-xks(:,i))*...
            (sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)+2*epsilon)/(sqrt(norm(x'-xks(:,i),2)^2+epsilon^2)+epsilon)^2;
    v = v + v_k;

end

function [u] = velocity_regularized(x,xks,fks)

% Same implementation as other velocity function, using regularized this
% time

u = zeros(2,1);
mu = .2;
n = size(xks,2);
epsilon = .3;

for i = 1:n
     uk = -fks(:,i)/(4*pi*mu) ...
         * ( log(sqrt(norm(x-xks(:,i)',2)^2 + epsilon^2)+2*epsilon) ...
            - epsilon * (sqrt(norm(x-xks(:,i)',2)^2 + epsilon^2)) ...
            /(sqrt(norm(x-xks(:,i)',2)^2 + epsilon^2)+epsilon ...
            *(sqrt(norm(x-xks(:,i)',2)^2 + epsilon^2)) ));
         + ...
         1/(4*pi*mu)*(dot(fks(:,i),x-xks(:,i)'))*(x-xks(:,i)) * ...
         (sqrt(norm(x-xks(:,i)',2)^2 + epsilon^2)+ 2*epsilon)/...
        ((sqrt(norm(x-xks(:,i)',2)^2 + epsilon^2)+epsilon)^2)*...
        sqrt(norm(x-xks(:,i)',2)^2);
     u = u + uk;
end

end
