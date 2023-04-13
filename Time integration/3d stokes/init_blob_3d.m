function [X] = init_blob(n)
% Init function returns location of [xi;yi] of n centers of each blob


delta = .01;
X = zeros(3,n);

for i = 1:n
    X(:,i) = [delta*(i-1); 0; 0];
end



