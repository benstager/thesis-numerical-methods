function [X] = Bp_3d(Xv,i) 
% returns ith blob position

X = [Xv(3*i-2);Xv(3*i-1);Xv(3*i)];
end