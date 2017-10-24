%% Time Integration

% This is a simple simulation of ODE solver
% Problem: y" = f --> y" + y = 0, y(0) = 2, y'(0) = 3
% Exact Solution: 2*cos(t) + 3*sin(t)
%
% Newmark time integration (Newmark, 1959) is used to integrate it   
%
% Copyright EYST <ezrayst@berkeley.edu>

clear
clc

dt = 0.1;
T = 8;
t = [0:dt:T];

%% Annotation: y is exact solution, u is the approximate solution
tt = [0:0.01:T];
y = 2*cos(tt) + 3*sin(tt);

%% ODE45 Solution
ff = @(t,y) [y(2) ; -y(1)];
[t_ode45, u_ode45] = ode45(ff, [0 T], [2; 3]);

%% Newmark Scheme
% Notation: u(:,1) = y, u(:,2) = y', u(:,3) = y"

u(1,1) = 2;
u(1,2) = 3;
u(1,3) = -u(1,1);
f(1) = -u(1,1);

beta = 0.25;
gamma = 0.5;

for ii = 1 : length(t) - 1;
    u(ii + 1, 1) = 1 / (1 + beta * dt^2) * ...
                   (u(ii, 1) + dt * u(ii, 2) + (1 - 2 * beta) / 2 * dt^2 * u(ii, 3));
    
    u(ii + 1, 2) = u(ii, 2) + (1 - gamma) * dt * u(ii, 3) - gamma * dt * u(ii + 1, 1);
    
    u(ii + 1, 3) = - u(ii + 1, 1);
end

%% Plotting
plot(tt,y,'k', t_ode45,u_ode45(:,1),'ob', t,u(:,1),'*r');
title('Solve y" = -y, y(0) = 2, y''(0) = 3'); 
xlabel('t');
ylabel('y');
legend('Exact', 'ODE45 Matlab', strcat('Newmark - \beta=',num2str(beta),', \gamma=',num2str(gamma), ', dt=', num2str(dt)), 'Location', 'SouthOutside');
