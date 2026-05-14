function [U, time_vec, Energy, snap_U, snap_t] = BCP2D_BHM_solver(dt, M1, tfinal, N, epsilon, alpha, snap_times)
%BCP2D_BHM_solver  2D BCP equation via BHM scheme
%   Returns U, energy history, and snapshots at snap_times

%% Domain
M_dom = 2;
a = 0; b = M_dom*pi;
h = (b-a)/N;

%% Grid and wavenumbers
x = gpuArray(a:h:b-h);

k = gpuArray([[0:N/2] [-N/2+1:-1]] ./ (M_dom/2));
[k1x, k1y] = meshgrid(k.^1, k.^1);
[kx, ky]   = meshgrid(k.^2, k.^2);
k2 = kx + ky;
k4 = k2.^2;
clear kx ky k

%% Initial condition
rng(1527, 'twister');
U = gpuArray(0.01*rand(N, N) + 0.25*0+.35);  %spheres or lamellar 

InitialMass = gather(sum(U(:)));

%% Parameters
eps2 = epsilon^2;
TheUbar = (1/(2*pi))^2 * h^2 * sum(sum(U));

%% LHS — Kim's scaling: CH + alpha implicit
lhs = 1 + dt*M1*k4*eps2 + dt*alpha;

%% Storage
nsteps   = round(tfinal/dt);
time_vec = zeros(1, nsteps+1);
Energy   = zeros(1, nsteps+1);
snap_U   = cell(length(snap_times), 1);
snap_t   = snap_times;

%% Initial energy
Energy(1)   = compute_energy_2D(U, k2, k1x, k1y, epsilon, alpha, TheUbar, h);
time_vec(1) = 0;

%% Time loop
hat_U = fft2(U);
t = 0; it = 0;

while t < tfinal - dt/2
    % Kim's scaling: CH + alpha*(ubar - u)
    RHS = eps2*(M1-1)*ifft2(k4.*fft2(U)) ...
        + ifft2(-k2.*fft2(U.^3 - U)) ...
        + alpha*TheUbar;
    hat_U = (hat_U + dt.*fft2(RHS)) ./ lhs;
    U = real(ifft2(hat_U));
    t = t + dt;
    it = it + 1;

    % Energy at every step
    Energy(it+1)   = compute_energy_2D(U, k2, k1x, k1y, epsilon, alpha, TheUbar, h);
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
fprintf('\nDone! t=%.1f  MassError=%.2e\n', t, abs(InitialMass-FinalMass)*h^2);

end

%% Energy function
function E = compute_energy_2D(U, k2, k1x, k1y, epsilon, alpha, Ubar, h)
    % Kim's BCP energy:
    % E = int F(phi) + (eps^2/2)|grad phi|^2 + (alpha/2)|grad psi|^2
    
    eps2 = epsilon^2;
    TheUbar = Ubar;
    
    % Double well
    dwell = (1/4)*U.^4 - (1/2)*U.^2;
    
    % Gradient energy: (eps^2/2)|grad u|^2
    dUx = real(ifft2(-1i*k1x.*fft2(U)));
    dUy = real(ifft2(-1i*k1y.*fft2(U)));
    grad_energy = (eps2/2)*(dUx.^2 + dUy.^2);
    
    % Nonlocal term: (alpha/2)|grad psi|^2
    psi_hat = fft2(U - TheUbar) ./ (k2 + (k2==0));
    psi_x = real(ifft2(-1i*k1x.*psi_hat));
    psi_y = real(ifft2(-1i*k1y.*psi_hat));
    nonlocal = (alpha/2)*(psi_x.^2 + psi_y.^2);
    
    E = gather(h^2 * sum(sum(dwell + grad_energy + nonlocal)));
end
