
lambda = .01;
Nts = 10.^(4:8);
y0 = 1;
tspan = [0,1];
dts = diff(tspan)./Nts;
error=zeros(1,length(Nts));
clock=zeros(1,length(Nts));

for i = 1:length(Nts)
    tic
    approx = eulerLin(lambda,tspan,y0,Nts(i),'dalquist',1,1);
    error(1,i) = abs(approx(end) - y0*exp(lambda*tspan(end)));
    clock(i) = toc;
end
loglog(dts,clock);
xlabel('time (sec)');
ylabel('error');
title('Error vs. Clock Time')