function U = heatExactSolution(xs,ts)
% returns exact solution to the heat equation at spatial points xs and time
% ts

ks = [1 2 10];
amps = [(8*pi/3) (4*pi/6) -(5*pi/6)];

U = zeros(size(xs));
for i = 1 : length(ks)
    U = U + amps(i) * sin(pi * ks(i) * xs) .* exp(-(pi * ks(i))^2 * ts);
end
end
