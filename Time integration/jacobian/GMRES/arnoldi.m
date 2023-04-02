function [h,q] = arnoldi(A,Q,i)

q = A*Q(:,i);

for j = 1:i
    h(j) = q' * Q(:,j);
    q = q - h(j) * Q(:,j);
end

h(i+1) = norm(q);
q = q / h(i+1);

end

