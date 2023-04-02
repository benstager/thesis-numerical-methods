
% using MATLAB functions to compute jacobian

syms x y z
A = jacobian([x*y*z,y^2,x + z],[x,y,z]);
J = zeros(size(A,1),size(A,2));
X = [1,2,3];

for i = 1:size(A,1)
    for j = 1:size(A,2)
        g = matlabFunction(A(i,j));
    end
end