function J = FD_JX(X,h)

% function returns the Jacobian given a system F and a vector X, stepsize h

F = {@(x,y) 2*x^2 - 3*y^2; @(x,y) 3*x^2 - 4*y^2};
J = zeros(size(F,1),length(X));

for i = 1:length(X)
    for j = 1:size(F,1)
        e = zeros(length(X),1);
        e(i) = 1;
        X_new = X + h*e;
        J(j,i) = (F{j}(X_new(1), X_new(2)) - F{j}(X(1),X(2)))/h;
    end
end

% redoing using generalized function notation

F = @(X) [2*X(1)^2-3*X(2)^2; 3*X(1)^2-4*X(2)^2];

for i = 1:length(X)
    e = zeros(length(X),1);
    e(i) = 1;
    J (:,i) = (F(X+h*e) - F(X))/h;
end

end

