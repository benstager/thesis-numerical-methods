
clear all;

% % system of equations: size i by 1 column vector
%
% f = {@(x,y) x^2+y^2; @(x,y) x^3-y^3};

syms x y
F = [2*x^2- 3*y^2; 3*x^2-4*y^2];
X = [x,y];
X0 = [1,2];
A = zeros(length(F), length(X));

for i = 1:length(F)
    for j = 1:length(X)
       g = matlabFunction(diff(F(i),X(j)));
       A(i,j) = g(X0(j));
    end
end
