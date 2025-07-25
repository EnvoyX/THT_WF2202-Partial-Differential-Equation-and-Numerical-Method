clear all; close all ; clc

% Discretization of Geometry
Length = 5;
Width = 5;

nelx = 10; % number of elements in x direction
nely = 10; % number of elements in y direction
nndx = nelx + 1; % number of nodes in x direction
nndy = nely + 1; % number of nodes in y direction

dx = linspace(0, Length, nndx);
dy = linspace(0, Width, nndy);

% Boundary Condition & Initial Condition
T = zeros(nndx, nndy); % Initialize solution for T

% Boundary Conditions from NIM 13123069
T_top = 131;
T_right = 23;
T_bottom = 6;
T_left = 9;


% Apply boundary conditions
T(1, :) = T_top;         % Top side
T(nndx, :) = T_bottom;      % Bottom side
T(:, 1) = T_left;        % Left side
T(:, nndy) = T_right;       % Right side

% Fix the corners by averaging overlapping BCs
T(1,1)   = (T_top + T_left) / 2;      % Top-left
T(1,nndy)   = (T_top + T_right) / 2;     % Top-right
T(nndx,1)   = (T_bottom + T_left) / 2;   % Bottom-left
T(nndx,nndy)   = (T_bottom + T_right) / 2;  % Bottom-right

% Initial Matrix T
disp(T);

% PLot initial and boundary conditions
figure;
[x,y] = meshgrid(dx,dy);
surf(x,y,T);
colorbar;
xlabel("X"); ylabel("Y"); zlabel("Temperature (°F)");
title('Initial and Boundary Conditions Temperature Distribution (Gauss-Seidel)');
view(45,30);

% Heatmap IC & BC  Plots
figure;
h = heatmap(T);
h.Title = 'Heated Plate Temperature Gradient Initial & Boundary Conditions';
h.XLabel = 'Left';
h.YLabel = 'Bottom'


% Numerical procedure
% Convergence criteria
epsilon_s = 1; % convergence tolerance in %
max_iter = 10000; % maximum iteration 

% Relaxition Paramter for faster calculation
omega = 1.5; % Relaxation Parameter

% Gauss-Seidel Algorithm
iter = 0;

while iter < max_iter
    T_prev_iter = T;
    for i = 2:nndx-1
        for j = 2:nndy-1
            T(i,j) = (1 - omega)*T(i,j) + omega*(0.25*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1)));
        end
    end

    % Convergence check
    if max(max(abs(T - T_prev_iter))) < epsilon_s
        fprintf('Converged in %d iterations with max error = %.4f%%\n', iter, max(max(abs(T - T_prev_iter))) );
        break
    end
    iter = iter + 1;
end



% Display final 3x3 interior values
fprintf('Final Temperature at Interior Points (%.f x %.f):\n',nelx-1,nely-1);
disp(T(2:nndx-1, 2:nndy-1));

% Optional: Full Grid Display
fprintf('Full Temperature Grid (%.f x %.f):\n',nelx,nely);
disp(T);

% Plotting the result
[X, Y] = meshgrid(dx, dy);
figure;
surf(X, Y, T, 'EdgeColor', 'none');
xlabel('X');
ylabel('Y');
zlabel('Temperature (°F)');
title('Steady-State Temperature Distribution (Gauss-Seidel)');
colorbar;
view(45,30);


% Heatmap Plots
figure;
h = heatmap(T);
h.Title = 'Heated Plate Temperature Gradient';
h.XLabel = 'Left';
h.YLabel = 'Bottom'

