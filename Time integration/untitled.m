

Nts = 16;
tspan = [0,1];
f = @(t,y) 3*y;
y0 = 1;

for i = 1:length(Nts)
    [ys,cpu] = euler(f,tspan,y0,Nts(i));
    ts = linspace(tspan(1),tspan(2),Nts(i)+1);
    plot(ts,ys, '-*', LineWidth = 6.0, color = 'black');
    xlabel('t'); ylabel('y(t)');
    hold on
end

hold on;
f = @(t) exp(3*t);
fplot(f,[0,1],LineWidth = 4.0, color = 'blue');
legend('Approximation','Exact');