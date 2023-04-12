function J = jac_FD(F,X,h)

J = zeros(length(X),length(X));

for i = 1:length(X)
    e = zeros(length(X),1);
    e(i) = 1;
    J(:,i) = (F(0,X+h*e) - F(0,X))/h;
end

end