% Input (known data)
L = 1.0;            % wall thickness in ft
alpha = 0.1;        % Diffusivity constant 
t_final = 1 ;             % Assume 1 Hour time
Ti = 100.0;         % initial temperature
Ts = 300.0;         % boundary temperature
dx = 0.05;          % Spatial step
dt = input('Enter time step(dt): ');  % Time Step
fprintf('Time step(dt): %d\n', dt);        
Nx = L/dx;          % number of spatial steps
r = alpha * dt / dx^2;

N = L/dx + 1 ; % Space nodes
M = t_final/dt + 1 ; % Time nodes

% Zero vectors of x & t
x = zeros(N, 1);
t = zeros(M, 1);


% Discretize space & time domain
for i = 1:N
    x(i) = 0 + (i - 1)*dx;  % step space domain
end

for n = 1:M
    t(n) = 0 + (n - 1)*dt;   % step time domain
end

% Outer Nodes
T_numerical = zeros(M,N);
T_numerical(:,1) = 300; %Left Boundary Condition
T_numerical(:,N) = 300; %Right Boundary Condition
T_numerical(1,2:N-1) = 100; %Initial Condition

% Fix the corners by averaging overlapping BCs
T_numerical(1,1)   = (300+100) / 2;   % Bottom-left
T_numerical(1,N)   = (300+100) / 2;  % Bottom-right


% Interior Nodes
for n= 1:M-1
    for i=2:N-1
        T_numerical(n+1, i) = r*T_numerical(n,i+1) + (1-2*r)*T_numerical(n,i) + r*T_numerical(n,i-1);
    end
end

% Plotting
x = linspace(0, L, N);
figure;
hold on;
for k = 1:size(T_numerical, 1)
    plot(x, T_numerical(k, :));
end
xlabel('Position (ft)');
ylabel('Temperature (F)');
title("Steady Heat Condunction in a Wall with time step (dt): " + num2str(dt) + "hr (FTCS Scheme)")
grid on;


% Comparison time (e.g., 1hr)
t_final = 1;
x = x';

% Compute analytical solution at t_final
n_terms = 100;                          % number of odd terms
n_vec = (1:2:2*n_terms)';               % odd integers: 1,3,5,...
C_n = -800 ./ (pi * n_vec);            % Coefficients
sin_terms = sin(n_vec * pi * x' / L);  % sin(n*pi*x)
exp_terms = exp(-alpha * (n_vec * pi / L).^2 * t_final);  % exp(-α(nπ/L)^2*t)
transient_sum = sum(C_n .* sin_terms .* exp_terms, 1);
T_analytical = Ts + transient_sum';    % Add 300 to get full solution


% Final Time Comparison Plot
n_idx = round(t_final / dt) + 1;  % Index in time for final step

figure;
plot(x, T_numerical(n_idx, :), 'r-', 'LineWidth', 3);
hold on;
plot(x, T_analytical, 'k--', 'LineWidth', 3);
xlabel('Position (ft)');
ylabel('Temperature (°F)');
title(sprintf('Final Time Comparison at t = %.2f hr', t_final));
legend('Numerical (FTCS)', 'Analytical');
grid on;
