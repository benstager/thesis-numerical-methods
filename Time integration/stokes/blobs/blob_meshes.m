function z = blob_meshes(x0,epsilon,xs,ys)

ny = length(ys);
nx = length(xs);
z = zeros(ny,nx);
blob = @(x,y) 3*epsilon^2/(2*pi*((norm([x;y]-x0,2)^2)+epsilon^2)^(5/2));

for i = 1:nx
    for j = 1:ny
        z(j,i) = blob(xs(i),ys(j));
    end
end


end

