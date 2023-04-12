function J = jac_ProdFD(F,X,V,h)

% computes matrix vector product Jv for X in R^j and small h
% of vector function system F
J(:) = (F(0,X+h*V) - F(0,X))/h;

end

