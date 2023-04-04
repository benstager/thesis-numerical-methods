function J = jac_ProdFD(F,X,V,h)

% computes matrix vector product Jv for X in R^j and small h
% of vector function system F
assert(length(X) == length(V), 'incorrect dimensions');

J = zeros(length(V),1);
J(:) = (F(X+h*V) - F(X))/h;

end

