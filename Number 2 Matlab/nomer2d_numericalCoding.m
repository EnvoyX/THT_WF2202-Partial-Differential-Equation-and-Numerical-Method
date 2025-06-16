% Parameters
alpha = 0.1;        % ft^2/hr
L = 1.0;            % wall thickness in ft
Ti = 100.0;         % initial temperature
Ts = 300.0;         % boundary temperature
dx = 0.05;
dt = 0.005;         % time step in hours
Nx = L/dx;            % number of spatial steps
r = alpha * dt / dx^2;

% Stability check
if r > 0.5
    error("r must be <= 0.5 for stability. r = %.3f", r);
end

% Time steps
Nt = 500;           % number of time steps

% Initialize temperature array
T = ones(1, Nx+1) * Ti;
T(1) = Ts;
T(end) = Ts;

% Store temperature for plotting every 50 steps
T_record = T;       % first time step

% Time-stepping loop
for n = 1:Nt
    T_new = T;
    for i = 2:Nx
        T_new(i) = T(i) + r * (T(i+1) - 2*T(i) + T(i-1));
    end
    T = T_new;
    T(1) = Ts;
    T(end) = Ts;

    if mod(n, 50) == 0
        T_record = [T_record; T];  % append new row
    end
end

% Plotting
x = linspace(0, L, Nx+1);
figure;
hold on;
for k = 1:size(T_record, 1)
    plot(x, T_record(k, :));
end
xlabel('Position (ft)');
ylabel('Temperature (F)');
title('Unsteady Heat Conduction in a Wall (FTCS)');
grid on;
