% Case 3
% Input (known data)
L = 1.0;            % wall thickness in ft
alpha = 0.1;        % Diffusivity constant 
T = 1 ;             % Assume 1 Hour time
Ti = 100.0;         % initial temperature
Ts = 300.0;         % boundary temperature
dx = 0.05;          % Spatial step
dt = 0.05;          % Time Step
Nx = L/dx;          % number of spatial steps
r = alpha * dt / dx^2;

N = L/dx + 1 ; % Space nodes
M = T/dt + 1 ; % Time nodes

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
U = zeros(M,N);
U(:,1) = 300; %Left Boundary Condition
U(:,N) = 300; %Right Boundary Condition
U(M,2:N-1) = 100; %Initial Condition
% Interior Nodes
for n= 1:M-2
    for i=2:N-1
        U(n+1, i) = r*U(n,i+1) + (1-2*r)*U(n,i) + r*U(n,i-1);
    end
end


disp(U)

% Plotting
x = linspace(0, L, Nx+1);
figure(1);
hold on;
for k = 1:size(U, 1)
    plot(x, U(k, :));
end
xlabel('Position (ft)');
ylabel('Temperature (F)');
title('Unsteady Heat Conduction in a Wall Case 3 (FTCS)');
grid on;


%fprintf('   t       x       U\n')
%for n= 1:M-2
    %for i=2:N-1
   %     U(n+1, i) = r*U(n,i+1) + (1-2*r)*U(n,i) + r*U(n,i-1);
  %      fprintf('%.4f   %.4f  %.4f\n' , t(n+1), x(i), U(n+1, i))
  %  end
 %   fprintf('   t       x       U\n')
%end

