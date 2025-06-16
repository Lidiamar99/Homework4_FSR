% Interpola i dati con passo regolare p.simTimeStep
t = (tout(1):p.simTimeStep:tout(end))';  % tempo uniforme

X  = interp1(tout, Xout,  t);
U  = interp1(tout, Uout,  t);
Xd = interp1(tout, Xdout, t);
Ud = interp1(tout, Udout, t);

% Costruzione del tempo per U2 a gradino
t2 = repelem(t, 2);
t2(1) = [];
t2(end+1) = t2(end);
U2 = repelem(U, 2, 1);

% Crea figura con 2 plot per riga (2x2)
fig = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact')

%% --- POSITION ---
ax1 = nexttile;
plot(t, X(:,1), 'r', t, X(:,2), 'g', t, X(:,3), 'b', ...
     t, Xd(:,1), 'r--', t, Xd(:,2), 'g--', t, Xd(:,3), 'b--', 'LineWidth', 1.2);
xlim([0 4.5]);
title('Position', 'FontSize', 20)
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Position [m]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12); grid on; box on;
legend({'$x$', '$y$', '$z$', '$x_{des}$', '$y_{des}$', '$z_{des}$'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'eastoutside');
expand_ylim(ax1)

%% --- VELOCITY ---
ax2 = nexttile;
plot(t, X(:,4), 'r', t, X(:,5), 'g', t, X(:,6), 'b', ...
     t, Xd(:,4), 'r--', t, Xd(:,5), 'g--', t, Xd(:,6), 'b--', 'LineWidth', 1.2);
xlim([0 4.5]);
title('Velocity', 'FontSize', 20)
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Velocity [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12); grid on; box on;
legend({'$v_x$', '$v_y$', '$v_z$', '$v_{x,des}$', '$v_{y,des}$', '$v_{z,des}$'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'eastoutside');
expand_ylim(ax2)

%% --- ANGULAR VELOCITY ---
ax3 = nexttile;
plot(t, X(:,16), 'r', t, X(:,17), 'g', t, X(:,18), 'b', ...
     t, Xd(:,16), 'r--', t, Xd(:,17), 'g--', t, Xd(:,18), 'b--', 'LineWidth', 1.2);
xlim([0 4.5]);
title('Angular Velocity', 'FontSize', 20)
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Angular Velocity [rad/s]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12); grid on; box on;
legend({'$\omega_x$', '$\omega_y$', '$\omega_z$', ...
        '$\omega_{x,des}$', '$\omega_{y,des}$', '$\omega_{z,des}$'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'eastoutside');
expand_ylim(ax3)

%% --- CONTROL FORCES (Fz) ---
ax4 = nexttile;
plot(t2, U2(:,3), 'r', t2, U2(:,6), 'g', t2, U2(:,9), 'b', t2, U2(:,12), 'k', ...
     t, Ud(:,3), 'r--', t, Ud(:,6), 'g--', t, Ud(:,9), 'b--', t, Ud(:,12), 'k--', 'LineWidth', 1.2);
xlim([0 4.5]);
title('Fz', 'FontSize', 20)
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Force [N]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 12); grid on; box on;
legend({'$U_2$ leg 1', '$U_2$ leg 2', '$U_2$ leg 3', '$U_2$ leg 4', ...
        '$U_d$ leg 1', '$U_d$ leg 2', '$U_d$ leg 3', '$U_d$ leg 4'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'eastoutside');
expand_ylim(ax4)

% Salvataggio automatico
if ~exist('plot_ex3', 'dir')
    mkdir('plot_ex3');
end
exportgraphics(fig, fullfile('plot_ex3', 'boundgait.pdf'));

% Funzione di espansione Y-axis
function expand_ylim(ax)
    yl = get(ax,'YLim');
    range = diff(yl);
    set(ax, 'YLim', yl + 0.05*range*[-1 1]);
end
