

Nts = [4 16 32];
tspan = [0,1];
f = @(t,y) 3*y;
y0 = 1;

for i = 1:length(Nts)
    [ys,cpu] = euler(f,tspan,y0,Nts(i));
    ts = linspace(tspan(1),tspan(2),Nts(i)+1);
    plot(ts,ys, '.-',MarkerSize = 35, LineWidth = 6.0);
    xlabel('t'); ylabel('y(t)');
    hold on
end

hold on;
f = @(t) exp(3*t);
fplot(f,tspan,LineWidth = 4.0, color = 'black');

legend('h = 1/4', 'h = 1/16', 'h = 1/32', 'Exact')


