
% Script plots the temporal error of centers for sequence of time steps
% taken to solve highly viscous stokes flow
clear all;

% grid parameters
n = 5;
X = init_blob(n);

% reference fine solution, typical of power length(Nts) + 1
Xv = reshape(X,[2*n,1]);
tspan = [0,1];
Nt_ref = 2^12;
[Xvs,cpu] = IMEXSemiLinearEulerGMRES(@f,tspan,Xv,Nt_ref);
Xvi_ref = Xvs(:,end);
% Xvi = reshape(Xvi,[2,length(Xvi)/2]);
% x_ref = Xvi(1,:);
% y_ref = Xvi(2,:);


methods = {@euler, @heun, @IMEXSemiLinearEulerGMRES};
Nts = 2.^(6:10);
dts = diff(tspan)./Nts;

error = zeros(length(Nts),length(methods));
time = zeros(length(Nts),length(methods));

% computing error at final time step for each cutoff position
for i = 1:length(methods)
    for j = 1:length(Nts)
        [Xvs,cpu] = methods{i}(@f,tspan,Xv,Nts(j));
        Xvi = Xvs(:,end);
%         Xvi = reshape(Xvi,[2,length(Xvi)/2]);
        x_coarse = Xvi(1,:);
        y_coarse = Xvi(2,:);
        % error(j,i) = norm([norm(x_ref-x_coarse,2);norm(y_ref-y_coarse,2)],2);
        error(j,i) = norm(Xvi_ref-Xvi);
        time(j,i) = cpu;
    end
end

getMethodName = @(f) functions(f).function;
legend_entries = cellfun(getMethodName, methods, 'UniformOutput', false);
set(0,'defaultAxesFontSize',13)

% convergence diagram
figure(1)
loglog(dts,error,'-*',LineWidth=4.0); hold on;
xlabel('step size (dt)');
ylabel('error');
title('Convergence Diagram for Stokes Flow');
legend(legend_entries); legend box on; 
ylim([10^-5, 10^2]);

% precision diagram
figure(2)
loglog(time,error,'-*',LineWidth=4.0);
xlabel('time (sec)');
ylabel('error');
title('Precision Diagram for Stokes Flow');
legend(legend_entries); legend box off;
ylim([10^-5, 10^2]);
        
