function U = heatExactSolution(xs,ts)
% returns exact solution to the heat equation at spatial points xs and time
% ts

U = (8*pi/3)*sin(pi*xs).*exp(-pi^2*ts);
end
