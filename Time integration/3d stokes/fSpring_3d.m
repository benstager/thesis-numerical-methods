function [F] = fSpring_3d(Xv,i,j)
% returns spring force on blob i, from j

k = 10^5; % spring constant
len_rest = .1;
len_curr = norm(Bp_3d(Xv,j)- Bp_3d(Xv,i),2);
f_mag = k*(len_curr-len_rest);
unit_vector = (Bp_3d(Xv,j) - Bp_3d(Xv,i))/len_curr;
F = f_mag * unit_vector;

end

