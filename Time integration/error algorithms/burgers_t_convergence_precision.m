% =========================================================================
% Burgers equation convergence and precision diagrams for time refinement.
% =========================================================================

clear;
close all;

% burgers

        Nx = 300; % spatial grid points
        Nt = 10.^(1:4); % timesteps to test
        Nt_ref = 10^5;
        eqn_name = sprintf('Burgers Equation (Nx = %i)', Nx);
        epsilon = .02;


        methods_classical = {
            @euler 
            @heun 
        };
        
        methods_IMEX = {
            @IMEXeuler
            @expIntegrator
            @IMEXmidpoint
            @IMEXdirk
            @expRK
        };
        
        Nc = length(methods_classical);
        Ni = length(methods_IMEX);

        [y0, xs, tspan] = burgersParameters(Nx);
        [f, L, N] = burgersOperators(Nx,epsilon);
        

        exact_IMEX  = burgersNumericalSolution(epsilon, Nt_ref, Nx, tspan);
%         approximate_IMEX = burgersNumericalSolution(epsilon, Nt_ref, 63, tspan);
%         spatialError = spatialError(approximate_IMEX(:,end), exact_IMEX(:,end));
        errorNorm = @(approx) norm(approx - exact_IMEX(:,end), 'inf');

%% run numerical experiments

error = zeros(length(Nt),Nc+Ni);
time = zeros(length(Nt),Nc+Ni);

for i = 1:Nc
    for j = 1:length(Nt)
        [ys,cpu_time] = methods_classical{i}(f,tspan,y0,Nt(j));
        error(j,i) = errorNorm(ys(:,end));
        time(j,i) = cpu_time;
    end
end

for i = 1:Ni
    for j = 1:length(Nt)
        [ys,cpu_time] = methods_IMEX{i}(L,N,tspan,y0,Nt(j));
        error(j,i+Nc) = errorNorm(ys(:,end));
        time(j,i+Nc) = cpu_time;
    end
end

dts = diff(tspan) ./ Nt;

%% generate plots 

% legend names for methods
methods_full = {
    @euler
    @heun
    @IMEXeuler
    @expIntegrator
    @IMEXmidpoint
    @IMEXdirk
    @expRK
    };

getMethodName = @(f) functions(f).function;
legend_entries = cellfun(getMethodName, methods_full, 'UniformOutput', false);
set(0,'defaultAxesFontSize',13)

% convergence diagram
figure(1)
loglog(dts,error,'-*',LineWidth=4.0);
%yline(spatialError, 'k--', LineWidth = 2.0);
ylim([1e-16, 1e2]);
xlabel('timestep (h)'); 
ylabel('error'); 
title(['Convergence Diagram - ', eqn_name]);
legend(legend_entries); legend box off;

% precision diagram
figure(2)
loglog(time,error,'-*',LineWidth=4.0);
% yline(spatialError, 'k--', LineWidth = 2.0);
ylim([1e-16,1e2])
xlabel('time (sec)'); 
ylabel('error'); 
title(['Precision Diagram - ', eqn_name]);
legend(legend_entries); legend box off;
