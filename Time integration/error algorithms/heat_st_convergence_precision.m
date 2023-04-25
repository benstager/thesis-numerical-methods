% =========================================================================
% Heat equation convergence and precision diagrams for space / time
% refinement.
% =========================================================================

clear;

% Test Spatial Grid Points 2^exps - 1
exps = 2:9;

% compute heat exact solution
Nx_ref = 2^(exps(end)) - 1;
[~,xs_ref,tspan] = heatParameters(Nx_ref);
exact_analytical = heatExactSolution(xs_ref(2:Nx_ref+1),tspan(end));

% methods{i} contains cell with { @integrator, @StepsizeSelector }
methods = {
    { @eulerLin,            @(n) 2 * diff(tspan) / ( (n + 1)^2 ) } % stability bound
    { @heunLin,             @(n) 2 * diff(tspan) / ( (n + 1)^2 ) } % stability bound
    { @backwardsEulerLin,   @(n) stepsizeSelector(tspan, n, 50, 4237, 2, 1) }   % obtained using Nx = 50 and [Userx,Usery] = ginput(1)
    { @impMidpointLin,      @(n) stepsizeSelector(tspan, n, 50, 80, 2, 2) }     % obtained using Nx = 50 and [Userx,Usery] = ginput(1)
    };

% FOR HEAT-MULTI-FREQ
methods = {
    { @eulerLin,            @(n) 2 * diff(tspan) / ( (n + 1)^2 ) } % stability bound
    { @heunLin,             @(n) 2 * diff(tspan) / ( (n + 1)^2 ) } % stability bound
    { @backwardsEulerLin,   @(n) stepsizeSelector(tspan, n, 50, 414, 2, 1) }   % obtained using Nx = 50 and [Userx,Usery] = ginput(1)
    { @impMidpointLin,      @(n) stepsizeSelector(tspan, n, 50, 18, 2, 2) }     % obtained using Nx = 50 and [Userx,Usery] = ginput(1)
    };

%% run numerical experiments

Nx = 2.^exps - 1;
error = zeros(length(Nx),length(methods));
time = zeros(length(Nx),length(methods));

for i = 1:length(methods)
    for j = 1:length(Nx)

        integrator = methods{i}{1};
        stepsizeSelector = methods{i}{2};

        Nt = ceil(diff(tspan)/stepsizeSelector(Nx(j)));
        [y0] = heatParameters(Nx(j));
        [A,f] = heatOperators(Nx(j));
        [ys,cpu_time] = integrator(A,tspan,y0,Nt);
        error(j,i) = spatialError(ys(:,end),exact_analytical);
        time(j,i) = cpu_time;

    end
end

%% Generate Plots

% legend names for methods
getMethodName = @(f) functions(f{1}).function;
legend_entries = cellfun(getMethodName, methods, 'UniformOutput', false);

% convergence diagram
figure(1)
loglog(Nx,error,LineWidth=3.0); hold on;
loglog(Nx, Nx.^(-2), 'k--'); hold off;
xlabel('grid size (Nx)');
ylabel('error');
title('Convergence Diagram for Heat Equation');
legend(legend_entries); legend box off;
ylim([10^-12,10^2]); xlim('tight');

% precision diagram
figure(2)
loglog(time,error,LineWidth=3.0);
xlabel('time (sec)');
ylabel('error');
title('Precision Diagram for Heat Equation');
legend(legend_entries); legend box off;
ylim([10^-16,10^2])