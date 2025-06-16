clc
close all
clear all

% Parameters
g = 9.81;          % gravity (m/s^2)
%%default values
%l = 1.00;           % leg length (m) %%MODIFY HERE%% 
%alpha = pi/8;      % half inter-leg angle (rad) %%MODIFY HERE%% 
%gamma = 0.08;       % slope angle (rad) %%MODIFY HERE%% 

%%only leg length variation
% l = 1.50;
% alpha=pi/8; 
% gamma = 0.08; 

% %only inter-leg angle variation
% l=1.00;
% alpha=pi/16;
% gamma=0.08;

% %only slope inclination variation
% l=1.00;
% alpha=pi/8;
% gamma=0.15;

%modified parameters in section 4.2.1
l = 0.9;
alpha = pi/12;
gamma = 0.08;

% Initial conditions
omega_1=sqrt(2*(g/l)*(1-cos(gamma-alpha)));
%thetadot0 must be greater than omega_1! the default thetadot0 was 0.95
%but it is less than omega_1 that is equal to 0.9754 
%thetadot0 = 0.98;
thetadot0 = 0.95;%default value
if (thetadot0 >= 0)
    theta0 = gamma-alpha;
else
    theta0 = gamma+alpha;
end

double_support = 0;

y0 = [theta0; thetadot0];

% Simulation settings
t0 = 0; %initial time
tf = 25; %final time
dt = 0.01; %max step time

% Time/state storage
T = [];
Y = [];

while t0 < tf
    options = odeset('Events', @(t, y) impact_event(t, y, alpha,gamma), 'MaxStep', dt);
    [t, y, te, ye, ie] = ode45(@(t, y) dynamics(t, y, g, l, double_support), [t0 tf], y0, options);
    
    T = [T; t];
    Y = [Y; y];

    if ~isempty(te)
        [y0,double_support] = impact_map(ye, alpha,g,l); % apply impact map
        t0 = te;
    else
        break;
    end
end

% Plot results
fig1 = figure('Renderer', 'painters', 'Position', [100 100 800 500]);
plot(T, Y(:,1), 'b', 'DisplayName', '$\theta$ [rad]', 'LineWidth', 1.5);
hold on;
plot(T, Y(:,2), 'r', 'DisplayName', '$\dot{\theta}$ [rad/s]', 'LineWidth', 1.5);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('States', 'Interpreter', 'latex', 'FontSize', 14);
title('Rimless Wheel Dynamics', 'FontSize', 18);
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
set(gca, 'FontSize', 12); grid on; box on;

% Save
exportgraphics(fig1, fullfile('plot_ex4', 'rimless_time_2nd2.pdf'));


fig2 = figure('Renderer', 'painters', 'Position', [100 100 800 500]);
plot(Y(:,1), Y(:,2), 'b', 'DisplayName', '$\theta$ vs $\dot{\theta}$', 'LineWidth', 1.5);
hold on;
plot(Y(1,1), Y(1,2), 'r*', 'LineWidth', 2.5, 'DisplayName', 'Initial Point');
xlabel('$\theta$ [rad]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 14);
title('Rimless Wheel Limit Cycle', 'FontSize', 18);
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
set(gca, 'FontSize', 12); grid on; box on;

% save
exportgraphics(fig2, fullfile('plot_ex4', 'rimless_phase_2nd2.pdf'));

function dydt = dynamics(~, y, g, l, ds)
    theta = y(1);
    thetadot = y(2);
    if (~ds)
        dtheta = thetadot;
        dthetadot = (g/l) * sin(theta);
    else
        dtheta = 0;
        dthetadot = 0;
    end
    dydt = [dtheta; dthetadot];
end

function [value, isterminal, direction] = impact_event(~, y, alpha,gamma)
    
    value = [y(1)-alpha-gamma; y(1)-gamma+alpha];% Trigger when theta = gamma+alpha
                                     %Trigger when theta = gamma-alpha
    isterminal = [1;1];         % Stop the integration
    direction = [1;-1];          % Detect only when increasing
end

function [yplus,ds] = impact_map(y_minus, alpha,g,l)%minus: before impact time; plus: after impact time
    if (y_minus(2)>=0)
        theta_plus = y_minus(1)-2*alpha;
    else
        theta_plus = y_minus(1)+2*alpha;
    end
    thetadot_plus = cos(2*alpha) * y_minus(2);
    if (thetadot_plus < 0.01*sqrt(g/l) && thetadot_plus >-0.01*sqrt(g/l)) 
        thetadot_plus = 0;
        ds = 1;
    else
        ds = 0;
    end
    yplus = [theta_plus; thetadot_plus];
end