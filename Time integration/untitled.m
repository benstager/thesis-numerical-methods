

Nts = [4 16];
tspan = [0,1];
f = @(t,y) -10*y;
y0 = 1;

for i = 1:length(Nts)
    [ys,cpu] = euler(f,tspan,y0,Nts(i));
    ts = linspace(tspan(1),tspan(2),Nts(i)+1);
    plot(ts,ys, '.-',MarkerSize = 35, LineWidth = 6.0);
    xlabel('t'); ylabel('y(t)');
    hold on
end

hold on;
f = @(t) exp(-15*t);
fplot(f,tspan,LineWidth = 4.0, color = 'black');
ylim([-2,2]);

legend('h = 1/4', 'h = 1/16', 'Exact');