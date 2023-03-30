function [X] = Bp(Xv,i) 
% returns ith blob position

X = [Xv(2*i-1);Xv(2*i)];
end