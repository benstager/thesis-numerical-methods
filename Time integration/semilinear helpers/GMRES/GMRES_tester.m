
% computing various solutions of Uxx x = b using GMRES

A = Uxx(1000,2/1000,true);
b = rand(1000,1);

% 1. tol 1e^-3

tol = 10^-3;
maxit = 1000;
x = gmres(A,b,[],tol,maxit);
plot(x);

% 2. tol 1e^-5

tol = 10^-5;
maxit = 1000;
x = gmres(A,b,[],tol,maxit);

% 3. Using a function handle 

tol = [10^-3 10^-5];
maxit = 1000;

for i = 1:length(tol)
    x = gmres(@(X) A*X,b,[],tol(i),maxit);
end





