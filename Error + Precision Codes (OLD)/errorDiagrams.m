tspan = [0,2];
y0 = 1;
A = -1;
B = 1;
epsilon = 1;
Nt = 10.^(3:8);
dts = diff(tspan)./Nt;
errorNorm = @(approx) abs(approx-exp(tspan(end)*A));
methods = {@eulerLin,@heunLin,@backwardsEulerLin,@impMidpointLin};

error = zeros(length(Nt),length(methods));
time = zeros(length(Nt),length(methods));

for j = 1:length(Nt)
    [ys,cpu_time] = eulerLin(A,B,tspan,y0,Nt(j),'dalquist',epsilon);
    error(j) = errorNorm(ys(:,end));
    time(j) = cpu_time;
end

loglog(time,error,LineWidth = 2.0);
xlabel('time (sec)')
ylabel('error')
title("Precision Diagram for Euler's Method")
legend('Euler')