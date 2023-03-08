close all;
clear all; %#ok<CLALL> 
set(gca,'Color','black')
p = nsidedpoly(1000, 'Center', [1 0], 'Radius', 1);
plot(p, 'FaceColor', 'red')
axis equal
hold on
plot([0 0], ylim,'black',LineWidth=1.5)               % Dashed Vertical Line at x=0
plot(xlim, [0 0],'black',LineWidth=1.5)    
title('Stability Region for Forward Euler')
xlabel('Re(z)')
ylabel('Im(z)') 