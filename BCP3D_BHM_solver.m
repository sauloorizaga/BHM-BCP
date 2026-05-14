function [U, time_vec, Energy, snap_U, snap_t] = BCP3D_BHM_solver(dt, M1, tfinal, N, epsilon, alpha, ubar, snap_times)
%BCP3D_BHM_solver  3D BCP equation via BHM scheme + GPU
%   snap_times: vector of times to save snapshots

%% Domain
M_dom = 2;
a = 0; b = M_dom*pi;
h = (b-a)/N;

%% Wavenumbers
k = gpuArray([[0:N/2] [-N/2+1:-1]] ./ (M_dom/2));
[k1x, k1y, k1z] = meshgrid(k.^1, k.^1, k.^1);
[kx,  ky,  kz]  = meshgrid(k.^2, k.^2, k.^2);
k2 = kx + ky + kz;
k4 = k2.^2;
clear kx ky kz k

%% Initial condition
rng(1527,'twister');
U = gpuArray(0.01*rand(N,N,N) + ubar);
InitialMass = gather(sum(U(:)));

%% Parameters
eps2 = epsilon^2;
TheUbar = (1/(2*pi))^3 * h^3 * sum(sum(sum(U)));

%% LHS — Kim's scaling
lhs = 1 + dt*M1*k4*eps2 + dt*alpha;

%% Storage
nsteps   = round(tfinal/dt);
time_vec = zeros(1, nsteps+1);
Energy   = zeros(1, nsteps+1);
snap_U   = cell(length(snap_times), 1);
snap_t   = snap_times;

%% Initial energy
Energy(1)   = compute_energy_3D(U, k2, k1x, k1y, k1z, epsilon, alpha, TheUbar, h);
time_vec(1) = 0;

%% Time loop
hat_U = fftn(U);
t = 0; it = 0;

while t < tfinal - dt/2
    RHS = eps2*(M1-1)*ifftn(k4.*fftn(U)) ...
        + ifftn(-k2.*fftn(U.^3 - U)) ...
        + alpha*TheUbar;
    hat_U = (hat_U + dt.*fftn(RHS)) ./ lhs;
    U = real(ifftn(hat_U));
    t = t + dt;
    it = it + 1;

    Energy(it+1)   = compute_energy_3D(U, k2, k1x, k1y, k1z, epsilon, alpha, TheUbar, h);
    time_vec(it+1) = t;

    % Save snapshots
    for s = 1:length(snap_times)
        if abs(t - snap_times(s)) < dt/2
            snap_U{s} = gather(U);
            fprintf('Snapshot saved at t=%.2f\n', t);
        end
    end
end

%% Gather
U = gather(U);
FinalMass = sum(U(:));
fprintf('\nDone! t=%.1f  MassError=%.2e\n', t, abs(InitialMass-FinalMass)*h^3);

end

%% Energy function — Kim's scaling
function E = compute_energy_3D(U, k2, k1x, k1y, k1z, epsilon, alpha, Ubar, h)
    eps2 = epsilon^2;

    % Double well
    dwell = (1/4)*U.^4 - (1/2)*U.^2;

    % Gradient energy
    dUx = real(ifftn(-1i*k1x.*fftn(U)));
    dUy = real(ifftn(-1i*k1y.*fftn(U)));
    dUz = real(ifftn(-1i*k1z.*fftn(U)));
    grad_energy = (eps2/2)*(dUx.^2 + dUy.^2 + dUz.^2);

    % Nonlocal
    psi_hat = fftn(U - Ubar) ./ (k2 + (k2==0));
    psi_x = real(ifftn(-1i*k1x.*psi_hat));
    psi_y = real(ifftn(-1i*k1y.*psi_hat));
    psi_z = real(ifftn(-1i*k1z.*psi_hat));
    nonlocal = (alpha/2)*(psi_x.^2 + psi_y.^2 + psi_z.^2);

    E = gather(h^3 * sum(sum(sum(dwell + grad_energy + nonlocal))));
end
