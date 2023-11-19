% iħ * Ψ_tt(r,t)= (-ħ² / 2m) * Ψ_rr(r,t) + V(r) * Ψ(r,t)
% http://www.niser.ac.in/~sbasak/p303_2010/14_15.09.pdf
clearvars -except f

% Parameters
L = 2e-8; % Width of the well in meters
N = 4000; % Number of spatial points
dr = L / (N - 1);
r = linspace(0, L, N)';

% Particle parameters
mass = 9.10938356e-31; % Mass of an electron in kg
hbar = 1.054571817e-34; % Reduced Planck's constant

q_e = 1.602176634e-19; % Elementary charge in Coulombs
epsilon0 = 8.854187817e-12; % Permittivity of free space

deltaE = (1/1^2 - 1/2^2) * (q_e^2 / (4 * pi * epsilon0))^2 * mass / (2 * hbar^2);
disp(deltaE);
disp(deltaE/q_e);
omega = deltaE / hbar; % Angular frequency of the transition
Tmn = 2 * pi / omega; % Period of the potential

% Time evolution parameters
dt = 2e-18; % Time step in seconds

% Time intervals in picoseconds
t1 = 0; %Tmn/4; % Time in the first excited state
t2 = 80*Tmn; % Time applying the force
t3 = 20*Tmn; % Time after the force is turned off
    
% Convert time intervals to steps
steps1 = floor(t1 / dt);
steps2 = floor(t2 / dt);
steps3 = floor(t3 / dt);

% Total number of time steps
num_steps = steps1 + steps2 + steps3;

% Finite difference representation of Laplacian (second derivative)
e = ones(N, 1);
Lap = spdiags([e -2*e e], -1:1, N, N) / dr^2;

% Kinetic energy operator (T = -hbar^2/(2m)*d^2/dr^2)
T = -hbar^2/(2*mass) * Lap;

% Initial wave function: ground state in an infinite well
a0 = 1.5*5.29177210903e-11; % Bohr radius
% Ground state
psi1 = 2*sqrt(1 / (pi * a0^3)) * exp(-r / a0) .* r;
% First excited state
psi2 = 1/4*sqrt(1 / (2 * pi * a0^3)) * (2 - r / a0) .* exp(-r / (2*a0)) .* r;
psi0=psi1;
psi0(1) = 0; psi0(end) = 0; % Enforce boundary conditions

% Preallocate psi over time
psi_t = zeros(N, num_steps);
psi_t(:, 1) = psi0;

V_r = -q_e^2/(4 * pi * epsilon0) ./ r;
V_r(1) = V_r(2);
V0=-V_r(1);

V_t = zeros(N, num_steps);
VV_t = zeros(N, num_steps);
VV_t(:, 1) = V0;
VV_t(1, 1) = V0 - V0 / 2;
VV_t(end, 1) = V0 + V0 /2;

V=0.25e-29;
disp(['t1 = ', num2str(steps1 * dt)]);
disp(['t2 = ', num2str((steps1 + steps2) * dt)]);
disp(['t3 = ', num2str((steps1 + steps2 + steps3) *dt)]);
% Time evolution
for t = 2:num_steps
    % During the second stage, apply the time-dependent potential
    if t > steps1 && t <= (steps1 + steps2)
        V_w = V*(cos(omega * (t * dt))) ./ r;
        V_w(1) = V_w(2);
        V_t(:, t) = V_r + V_w;
    else
        V_t(:, t) = V_r;
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
hPotential = plot(r, VV_t(:, 1), 'LineWidth', 2);
xlabel('Position (m)');
ylabel('V(r)');
title('Potential Energy at t=0');
axis([0 L -V0/100 V0/100]);
grid on;

% Subplot for the modulus of the wave function
subplot(4, 1, 2);
hModulus = plot(r, abs(psi_t(:, 1)), 'LineWidth', 2);
hold on; % Keep the current plot
hPsi2 = plot(r, abs(2.5*psi2(:, 1)), 'LineWidth', 2); % Plotting psi2
hold off; % Release the plot for next plots

xlabel('Position (m)');
ylabel('|Ψ|');
title('Modulus of Wave Function at t=0');
axis([0 L/10 0 max(abs(psi_t(:)))]);
grid on;

% Subplot for the phase of the wave function
subplot(4, 1, 3);
hPhase = plot(r, angle(psi_t(:, 1)), 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Phase (radians)');
title('Phase of Wave Function at t=0');
axis([0 L/10 -pi pi]);
grid on;

subplot(4, 1, 4);
hProbDensity = plot(r, abs(psi_t(:, 1)).^2, 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Probability Density');
title('Probability Density at t=0');
axis([0 L/10 0 max(abs(psi_t(:)).^2)]); % Changed to max of entire psi_t for proper scaling
grid on;

% Animation loop
for idx = 2:10:num_steps
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
