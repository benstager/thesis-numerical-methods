function J = JF(X)

% function returns Jacobian matrix at point X in R^j

J_fun = {@(x,y) 4*x, @(x,y) -6*y; @(x,y) 6*x, @(x,y) -8*y};
J = zeros(size(J_fun,1),length(X));

assert(length(X) == size(J_fun,2), 'x_j does not equal input vector');

for i = 1:length(X)
    for j = 1:size(J_fun,1)
        J(j,i) = J_fun{j,i}(X(1),X(2));
    end
end


