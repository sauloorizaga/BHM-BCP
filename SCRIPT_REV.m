%% BCP 3D - Main Script  
tic;

%% USER PARAMETERS
N       = 128/2;
dt      = 0.01;
M1      = 7.5;
tfinal  = 100;
epsilon = 0.1;
alpha   = 15;
%ubar    = 0.15;     % 
snap_times = [5, 10, 20, 100];   ubar    = 0.0;   %Lamellar
snap_times = [20, 25, 30, 100];   ubar    = 0.17;   %Gyroid

%% RUN SOLVER
[U, time_vec, Energy, snap_U, snap_t] = BCP3D_BHM_solver(dt, M1, tfinal, N, epsilon, alpha, ubar, snap_times);

%% DOMAIN
a = 0; b = 2*pi; h = (b-a)/N;
x = a:h:b-h;
[X,Y,Z] = meshgrid(x,x,x);

snap_labels = {['$t=' num2str(snap_times(1)) '$'], ...
               ['$t=' num2str(snap_times(2)) '$'], ...
               ['$t=' num2str(snap_times(3)) '$'], ...
               ['$t=' num2str(snap_times(4)) '$']};

%% FIGURES 1-4: Individual isosurfaces
for s = 1:4
    figure(s); clf;
    isosurface(X,Y,Z,snap_U{s},-.15);
    isosurface(X,Y,Z,snap_U{s},-.05);
    isosurface(X,Y,Z,snap_U{s},.05);
    isosurface(X,Y,Z,snap_U{s},.15);
    ax = gca; ax.FontSize = 14;
    camlight; lighting phong;
    axis([a b a b a b]);
    title(snap_labels{s},'Interpreter','latex','FontSize',18);
    saveas(gcf, ['BCP3D_snap' num2str(s) '_ubar' num2str(ubar*100) '.png']);
end

%% FIGURE 5: Energy — clean, no dots, no insets
figure(5); clf;
%set(gcf, 'Position', [100 100 1400 500]);
semilogx(time_vec, Energy, 'b-', 'LineWidth', 2.5);
grid on; grid minor;
ax = gca; ax.FontSize = 14;
ax.TickLabelInterpreter = 'latex';
xlabel('Time ($t$)','Interpreter','latex','FontSize',16);
ylabel('Free Energy $\mathcal{E}(u)$','Interpreter','latex','FontSize',16);
title(['BCP 3D Energy -- BHM Scheme, $\bar{u}=' num2str(ubar) '$'], ...
    'Interpreter','latex','FontSize',18);
axis tight;
saveas(gcf, ['BCP3D_energy_ubar' num2str(ubar*100) '.png']);

a_time = toc;
minutes_hours = [a_time/60, a_time/3600]