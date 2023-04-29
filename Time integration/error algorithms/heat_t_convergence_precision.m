% =========================================================================
% Heat equation convergence and precision diagrams for time refinement.
% =========================================================================

clear;
close all;

eq = 'heat'; % 'heat' or 'dalquist'

switch eq
    
    case 'heat'
        
        Nx = 200; % spatial grid points
        Nt = 10.^(2:6); % timesteps to test
        eqn_name = sprintf('Heat Equation (Nx = %i)', Nx);
        
        methods = {
            @eulerLin 
            @heunLin 
            @backwardsEulerLin
            @impMidpointLin
        };
    
        [y0, xs, tspan] = heatParameters(Nx);
        A = heatOperators(Nx);
    
        exact_analytical = heatExactSolution(xs(2:Nx+1),tspan(end));
        exact_numerical  = expm(tspan(end)*A)*y0;
  
        spatialError = norm(exact_analytical - exact_numerical, 'inf');
        errorNorm = @(approx) norm(approx - exact_numerical, 'inf');
   
    case 'dalquist' % y' = A y
        
        A = 3;        
        Nt = 10.^(3:7);
        tspan = [0,5];
        y0 = 1;
        eqn_name = sprintf('Dalquist');
        
        exact_analytical = exp(tspan(end)*A) * y0;
        
        methods = {@eulerLin
                   @heunLin
                   @backwardsEulerLin
                   @impMidpointLin
            };
        errorNorm = @(approx) abs(approx - exact_analytical);
        spatialError = 0;
       
end

%% run numerical experiments

error = zeros(length(Nt),length(methods));
time = zeros(length(Nt),length(methods));

for i = 1:length(methods)
    for j = 1:length(Nt)
        [ys,cpu_time] = methods{i}(A,tspan,y0,Nt(j));
        error(j,i) = errorNorm(ys(:,end));
        time(j,i) = cpu_time;
    end
end

dts = diff(tspan) ./ Nt;

%% generate plots 

% legend names for methods
getMethodName = @(f) functions(f).function;
legend_entries = cellfun(getMethodName, methods, 'UniformOutput', false);
set(0,'defaultAxesFontSize',13)

% convergence diagram
figure(1)
loglog(dts,error,LineWidth=4.0);
%yline(spatialError, 'k--', LineWidth = 2.0);
ylim([1e-16, 1e2]);
xlabel('timestep (h)'); 
ylabel('error'); 
title(['Convergence Diagram - ', eqn_name]);
legend(legend_entries); legend box off;

% precision diagram
figure(2)
loglog(time,error,LineWidth=4.0);
%yline(spatialError, 'k--', LineWidth = 2.0);
ylim([1e-16,1e2])
xlabel('time (sec)'); 
ylabel('error'); 
title(['Precision Diagram - ', eqn_name]);
legend(legend_entries); legend box off;
