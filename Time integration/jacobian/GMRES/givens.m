function [h,cs_i,sn_i] = givens(h,cs,sn,i)

for j = 1:i
    temp   =  cs(j) * h(j) + sn(j) * h(j + 1);
    h(j+1) = -sn(j) * h(j) + cs(j) * h(j + 1);
    h(j)   = temp;
end

[cs_i, sn_i] = givens_set(h(i), h(i + 1));

h(i) = cs_i * h(i) + sn_i * h(i + 1);
h(i + 1) = 0.0;


end

