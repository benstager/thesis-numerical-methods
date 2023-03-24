function v = velocity(x,xks,fks)

v = zeros(2,1);
mu = 1;
n = size(xks,2);


for i = 1:n
    v_k = (-fks(:,i)/(4*pi*mu))*log(norm(x'-xks(:,i),2))...
    +(dot(fks(:,i),x'-xks(:,i)))*(x'-xks(:,i))/(4*pi*mu*norm(x'-xks(:,i),2)^2);
    v = v + v_k;
end

