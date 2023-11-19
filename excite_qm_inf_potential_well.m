% iħ * Ψ_tt(x,t)= (-ħ² / 2m) * Ψ_xx(x,t) + V(x) * Ψ(x,t)

% Parameters
L = 1e-8; % Width of the well in meters
N = 5000; % Number of spatial points
dx = L / (N - 1);
x = linspace(0, L, N)';

% Particle parameters
mass = 9.10938356e-31; % Mass of an electron in kg
hbar = 1.054571817e-34; % Reduced Planck's constant


deltaE = (9 - 1) * pi^2 * hbar^2 / (2 * mass * L^2); % Energy difference between n=2 and n=1
omega = deltaE / hbar; % Angular frequency of the transition
T = 2 * pi / omega; % Period of the potential

% A reasonable starting guess for V0 could be a fraction of the energy difference
V0 = 0.25 * deltaE; % This is just a guess; needs to be fine-tuned

% Time evolution parameters
dt = 1e-16; % Time step in seconds


% Time intervals in picoseconds
t1 = T/4; % Time in the first excited state
t2 = 12*T; % Time applying the force
t3 = 14*T; % Time after the force is turned off
    
% Convert time intervals to steps
steps1 = floor(t1 / dt);
steps2 = floor(t2 / dt);
steps3 = floor(t3 / dt);

% Total number of time steps
num_steps = steps1 + steps2 + steps3;

% Finite difference representation of Laplacian (second derivative)
e = ones(N, 1);
Lap = spdiags([e -2*e e], -1:1, N, N) / dx^2;

% Kinetic energy operator (T = -hbar^2/(2m)*d^2/dx^2)
T = -hbar^2/(2*mass) * Lap;

% Define the perturbation potential parameters
% V0 = 0.5e-20; % Amplitude of the potential
deltaE = (pi^2 * hbar^2 / (2 * mass * L^2)) * (9 - 1); % Energy diff between second and first excited state
omega = deltaE / hbar; % Resonant angular frequency

% Initial wave function: ground state in an infinite well
psi0 = sqrt(2/L) * sin(pi*x/L);
psi0(1) = 0; psi0(end) = 0; % Enforce boundary conditions

% Preallocate psi over time
psi_t = zeros(N, num_steps);
psi_t(:, 1) = psi0;
V_t = zeros(N, num_steps);
VV_t = zeros(N, num_steps);
VV_t(:, 1) = V0;
VV_t(1, 1) = V0 - V0 / 2;
VV_t(end, 1) = V0 + V0 /2;
disp(['t1 = ', num2str(steps1 * dt)]);
disp(['t2 = ', num2str((steps1 + steps2) * dt)]);
disp(['t3 = ', num2str((steps1 + steps2 + steps3) *dt)]);
% Time evolution
for t = 2:num_steps
    % During the second stage, apply the time-dependent potential
    if t > steps1 && t <= (steps1 + steps2)
          V_t(:, t) = V0 * sin(pi * x / L) .* cos(omega * (t * dt));
    else
        V_t(:, t) = zeros(N, 1);
    end
    VV_t(:, t) = V_t(:, t);
    VV_t(1, t) = 0;
    VV_t(end, t) = 0;
    V_t(1, t) = Inf; 
    V_t(end, t) = Inf; % Infinite potential at the boundaries

    % Update the Hamiltonian with the time-dependent potential
    H_t = T + spdiags(V_t(:,t), 0, N, N);

    % Crank-Nicolson matrices with the time-dependent potential
    A = speye(N) - 1i * dt / (2 * hbar) * H_t;
    B = speye(N) + 1i * dt / (2 * hbar) * H_t;
    A = A(2:end-1, 2:end-1); % Remove the infinities
    B = B(2:end-1, 2:end-1);

    % Update the wave function using Crank-Nicolson with the potential
    psi_interior = A \ (B * psi_t(2:end-1, t-1));
    psi_t(:, t) = [0; psi_interior; 0]; % Enforce boundary conditions
end

% Check if 'f' exists and is a valid figure handle, otherwise create a new figure
if ~exist('f', 'var') || ~isvalid(f)
    f = figure;
else
    figure(f); % Brings the existing figure with handle 'f' to the front
end

% Subplot for the potential energy
subplot(4, 1, 1);
hPotential = plot(x, VV_t(:, 1), 'LineWidth', 2);
xlabel('Position (m)');
ylabel('V(x)');
title('Potential Energy at t=0');
axis([0 L -V0 V0]);
grid on;

% Subplot for the modulus of the wave function
subplot(4, 1, 2);
hModulus = plot(x, abs(psi_t(:, 1)), 'LineWidth', 2);
xlabel('Position (m)');
ylabel('|Ψ|');
title('Modulus of Wave Function at t=0');
axis([0 L 0 max(abs(psi_t(:)))]);
grid on;

% Subplot for the phase of the wave function
subplot(4, 1, 3);
hPhase = plot(x, angle(psi_t(:, 1)), 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Phase (radians)');
title('Phase of Wave Function at t=0');
axis([0 L -pi pi]);
grid on;

subplot(4, 1, 4);
hProbDensity = plot(x, abs(psi_t(:, 1)).^2, 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Probability Density');
title('Probability Density at t=0');
axis([0 L 0 max(abs(psi_t(:)).^2)]); % Changed to max of entire psi_t for proper scaling
grid on;

% Animation loop
for idx = 2:20:num_steps
    % Update the potential energy plot
    set(hPotential, 'YData', VV_t(:, idx));
    % Update the wave function title with current time
    title(subplot(4,1,1), ['Potential Energy at t=', num2str(idx*dt), ' ns']);

    % Update the wave function module plot
    set(hModulus, 'YData', abs(psi_t(:, idx)));
    % Update the wave function title with current time
    title(subplot(4,1,2), ['Wave Function at t=', num2str(idx*dt), ' ns']);

    % Update the wave function phase plot
    set(hPhase, 'YData', angle(psi_t(:, idx)));
    % Update the wave function title with current time
    title(subplot(4,1,3), ['Wave Function at t=', num2str(idx*dt), ' ns']);

    % Update the probability density plot
    set(hProbDensity, 'YData', abs(psi_t(:, idx)).^2);
    % Update the probability density title with current time
    title(subplot(4,1,4), ['Probability Density at t=', num2str(idx*dt), ' ns']);

    drawnow; % Refresh the plot

    % Pause for a brief moment to control the speed of the animation
    pause(0.01);

end
