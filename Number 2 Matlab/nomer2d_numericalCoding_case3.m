% Case 3
% Input (known data)
L = 1.0;            % wall thickness in ft
alpha = 0.1;        % Diffusivity constant 
time = 1 ;          % Assume 1 Hour time
Ti = 100.0;         % initial temperature
Ts = 300.0;         % boundary temperature
dx = 0.05;          % Spatial step
dt = 0.05;          % Time Step
Nx = L/dx;          % number of spatial steps
r = alpha * dt / dx^2;

N = L/dx + 1 ; % Space nodes
M = time/dt + 1 ; % Time nodes

% Zero vectors of x & t
x = zeros(N, 1);
t = zeros(M, 1);


% Discretize space & time domain
for i = 1:N
    x(i) = 0 + (i - 1)*dx;
end

for n = 1:M
    t(n) = 0 + (n - 1)*dt;
end

% Outer Nodes
T = zeros(M,N);
T(:,1) = 300; %Left Boundary Condition
T(:,N) = 300; %Right Boundary Condition
T(1,2:N-1) = 100; %Initial Condition

% Fix the corners by averaging overlapping BCs
T(1,1)   = (300+100) / 2;   % Bottom-left
T(1,N)   = (300+100) / 2;  % Bottom-right

% Interior Nodes
for n= 1:M-1
    for i=2:N-1
        T(n+1, i) = r*T(n,i+1) + (1-2*r)*T(n,i) + r*T(n,i-1);
    end
end
% Result of Matrix T
disp(T)

% Plotting
x = linspace(0, L, N);
figure(1);
hold on;
for k = 1:size(T, 1)
    plot(x, T(k, :));
end
xlabel('Position (ft)');
ylabel('Temperature (F)');
title('Unsteady Heat Conduction in a Wall Case 3 (FTCS)');
grid on;
