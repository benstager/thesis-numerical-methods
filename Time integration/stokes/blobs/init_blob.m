function [X] = init_blob(n)
% blob function returns location of [xi,yi] of n centers of each blob

epsilon = .1;
delta = .17;
X = zeros(2,n);


for i = 1:n
    X(:,i) = [delta*(i-1); 0];
end



