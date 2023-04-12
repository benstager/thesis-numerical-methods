function x = GMRES_scratch(A,b,x,m,epsilon)
% function GMRES_scratch attempts to return solution x to Ax = b

% A - square matrix
% b - solution vector
% x - initial guess and resulting solution
% m - iterations
% epsilon - threshold for b - Ax < epsilon

n = length(A);
b_norm = norm(b);
r = b - A*x;
Q(:,1) = (1/norm(r))*r;
e1 = zeros(m+1,1);
e1(1) = 1;
beta = norm(r) * e1;
sn = zeros(m, 1);
cs = zeros(m, 1);

for i = 1:m

    [H(1:i+1,i),Q(:,i+1)] = arnoldi(A,Q,i);
    [H(1:i+1,i),cs(i),sn(i)] = givens(H(1:i+1,i),cs,sn,i);
    
    beta(i+1) = -sn(i) * beta(i);
    beta(i) = cs(i) * beta(i);

end

y = H(1:i, 1:i) \ beta(1:i);
x = x + Q(:, 1:i) * y;

end


