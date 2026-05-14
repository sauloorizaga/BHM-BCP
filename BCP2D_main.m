%% BCP 2D - Main Script
tic;
 
%% USER PARAMETERS
N       = 128;
dt      = 0.01;
M1      = 10;
tfinal  = 100;
epsilon = 0.1/2;
alpha   = 15;
 
snap_times = [0.1, 1, 10, 100];   % snapshot times
 
%% RUN SOLVER
[U, time_vec, Energy, snap_U, snap_t] = BCP2D_BHM_solver(dt, M1, tfinal, N, epsilon, alpha, snap_times);
 
%% PLOTTING
a = 0; b = 2*pi; h = (b-a)/N;
x = linspace(a, b-h, N);
[X, Y] = meshgrid(x, x);
 
%--- Figure 1: 4 snapshots panel ---
figure(1); clf;
set(gcf, 'Position', [100 100 1200 350]);
title(titles{s}, 'Interpreter','latex','FontSize',18);
for s = 1:4
    subplot(1,4,s);
    pcolor(X, Y, snap_U{s}); shading interp;
    colormap(parula); 
    axis equal tight off;
    title(titles{s}, 'Interpreter','latex','FontSize',18);
end
sgtitle('BCP 2D Morphology Evolution -- BHM Scheme', 'Interpreter','latex','FontSize',18);
 
%--- Figure 2: Energy curve with embedded snapshot markers ---
figure(2); clf;
semilogx(time_vec, Energy, 'b-', 'LineWidth', 2.5);
hold on;
grid on; grid minor;
ax = gca; ax.FontSize = 18;
ax.TickLabelInterpreter = 'latex';
xlabel('Time ($t$)',                      'Interpreter','latex','FontSize',18);
ylabel('Free Energy $\mathcal{E}(u)$',    'Interpreter','latex','FontSize',18);
title('BCP 2D Energy Dissipation -- BHM Scheme', 'Interpreter','latex','FontSize',18);
axis tight;
 
% Mark snapshot times on energy curve
for s = 1:length(snap_times)
    [~, idx] = min(abs(time_vec - snap_times(s)));
    plot(time_vec(idx), Energy(idx), 'ko', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'k', 'LineWidth', 2);
    text(time_vec(idx), Energy(idx)*1.15, ['$t=' num2str(snap_times(s)) '$'], ...
        'Interpreter','latex','FontSize',12,'HorizontalAlignment','center');
end
 
%--- Figure 3: Final pcolor ---
figure(3); clf;
pcolor(X, Y, U); shading interp;
colormap(parula); colorbar;
ax = gca; ax.FontSize = 18;
axis equal tight;
title(['BCP 2D, $N=' num2str(N) '$, $t=' num2str(tfinal) '$, $\alpha=' num2str(alpha) '$'], ...
    'Interpreter','latex','FontSize',18);

%insertin figures in E plot
%--- Figure 10: Energy + embedded snapshot insets ---
figure(10); clf;
set(gcf, 'Position', [100 100 1400 500]);   % wide figure
semilogx(time_vec, Energy, 'b-', 'LineWidth', 2.5);
hold on;
grid on; grid minor;
ax = gca; ax.FontSize = 18;
ax.TickLabelInterpreter = 'latex';
xlabel('Time ($t$)', 'Interpreter','latex','FontSize',18);
ylabel('Free Energy $\mathcal{E}(u)$','Interpreter','latex','FontSize',18);
title('BCP 2D Energy Dissipation -- BHM Scheme','Interpreter','latex','FontSize',18);
axis tight;

% Mark 4 circles on curve
for s = 1:length(snap_times)
    [~, idx] = min(abs(time_vec - snap_times(s)));
    plot(time_vec(idx), Energy(idx), 'ko', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'k', 'LineWidth', 2);
end

% 4 inset axes — positioned roughly, drag them in figure editor! Labyrinth
inset_pos = [0.16 0.65 0.12 0.25;   % t=0.01
             0.32 0.18 0.12 0.25;   % t=1 — bit higher (was 0.10)
             0.55 0.32 0.12 0.25;   % t=10
             0.77 0.18 0.12 0.25];  % t=100
         
         
inset_pos = [0.16 0.65 0.12 0.25;   % t=0.01       SPHERES
             0.37 0.45 0.12 0.25;   % t=1 — higher
             0.58 0.28 0.12 0.25;   % t=10 — lower and right
             0.77 0.18 0.12 0.25];  % t=100
         

for s = 1:4
    ax_in = axes('Position', inset_pos(s,:));
    pcolor(X, Y, snap_U{s}); shading interp;
    colormap(parula); axis equal tight off;
end
 
a_time = toc;
minutes_hours = [a_time/60, a_time/3600]