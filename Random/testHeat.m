
N = 100;

xs = linspace(-1,1,N+2)';
dx = xs(2) - xs(1);
xs([1,end]) = [];

f = @(x) sin(pi*x);
fpp = @(x) -pi^2*sin(pi*x);

fv = f(xs);
fvpp = fpp(xs);

Dxx = Uxx(N,dx);
fvpp_approx = Dxx*fv;


subplot(2,1,1);
plot(xs,fvpp);
hold on
plot(xs,fvpp_approx);

subplot(2,1,2)
plot(fvpp_approx - fvpp);
%%
k = 20;
Ns = 2.^(1:15);
dxs = 2./(Ns+1);
f = @(x) sin(pi*k*x);
fpp = @(x) -pi^2*k^2*sin(k*pi*x);
f = @(x) 5*exp(-5*x^2);
error = zeros(1,length(Ns));

for i = 1:length(Ns)
    xs = linspace(-1,1, Ns(i)+2)';
    xs([1,end]) = [];


    Dxx = Uxx(Ns(i),dxs(i));
    fvpp = fpp(xs);
    fv = f(xs);
    fvpp_approx = Dxx*fv;
    error(i) = norm(fvpp_approx - fvpp,inf);
end



loglog(dxs,error,'--.');
hold on
loglog(dxs,dxs.^2,'red');
