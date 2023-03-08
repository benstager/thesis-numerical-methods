% Euler: usual formula
% Heun: usual formula 
% IMEX: referenceBurgers(epsilon,10^4,127,tspan), Nt0 = 100, Nx0 = 63
% ExpInt: referenceBurgers(epsilon,10^4,127,tspan), Nt0 = 100, Nx0 = 63

% # of interior spatial points
Ns = 2.^(4:7) - 1;
Ns_ref = 2^8-1;
Nt_ref = 5000;
epsilon = .02;
% amount of temporal points based on spatial grid
methods = {@euler,@heun};
methodsIMEX = {@IMEXeuler, @expIntegrator};

% What is p for IMEX and Exp? Should it be 2?
% epsilon here or no?
dts = {@(n) 2/(epsilon*(n+1)^2), @(n) 2/(epsilon*(n+1)^2),@(n) diff(tspan)/((100)/(63/n)^2),@(n) diff(tspan)/((100)/(63/n)^2)};
tspan = [0,.4];
xspan = [-1,1];
dxs = diff(xspan)./Ns;

% Should method be dynamic or static here?
exact = referenceBurgers(epsilon,Nt_ref,Ns_ref,tspan);
Nm = length(methods);
Nmi = length(methodsIMEX);
time = zeros(length(Ns),Nm + Nmi);
error = zeros(length(Ns),Nm+Nmi);


for i = 1:Nm
    for j = 1:length(Ns)
        Nt = ceil(diff(tspan)/dts{i}(Ns(j)));
        y0 = burgersIC(Ns(j));
        f = burgersOperators(epsilon,Ns(j));
        [ys,cpu_time] = methods{i}(f,tspan,y0,Nt);
        error(j,i) = spatialError(ys(:,end),exact(:,end));
        time(j,i) = cpu_time;
    end
end

for i= 1:Nmi
    for j = 1:length(Ns)
        Nt = ceil(diff(tspan)/dts{i}(Ns(j)));
        y0 = burgersIC(Ns(j));
        [f,L,Nl] = burgersOperators(epsilon,Ns(j));
        [ys,cpu_time] = methodsIMEX{i}(L,Nl,tspan,y0,Nt);
        error(j,i+2) = spatialError(ys(:,end),exact(:,end));
        time(j,i+2) = cpu_time;
    end
end


methods = {@euler,@heun,@IMEXeuler,@expIntegrator};
getMethodName = @(f) functions(f).function;
legend_entries = cellfun(getMethodName, methods, 'UniformOutput', false);

figure(1)
loglog(Ns,error,LineWidth=2.0);
xlabel('Grid size (Nx)'); ylabel('error'); title('test');
legend(legend_entries);
ylim([10^-12,10^2])

figure(2)
loglog(time,error,LineWidth=2.0);
xlabel('time (sec)'); ylabel('error'); title('Precision Diagram for Burgers Equation);
legend(legend_entries);
ylim([10^-12,10^2])