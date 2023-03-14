function [fns] = forces(xn)
%FORCES Summary of this function goes here
%   Detailed explanation goes here
fns = zeros(size(xn));
fns(1,:) = 1;
fns(2,:) = -1;

end

