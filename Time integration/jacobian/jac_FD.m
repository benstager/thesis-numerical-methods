function J = jac_FD(F,X,h)

% computes numerical jacobian at point X in R^j for small h
% for function F, with f_1...f_i

for i = 1:length(X)
    e = zeros(length(X),1);
    e(i) = 1;
    J(:,i) = (F(X+h*e) - F(X))/h;
end

end

