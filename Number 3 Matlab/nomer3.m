clc; clear all; close all;

% Grid size (number of interior nodes)
N = 5;
T = zeros(N+2, N+2);  % Include boundaries

% Boundary conditions from NIM 13123069
T_top = 131;
T_right = 23;
T_bottom = 6;
T_left = 9;

% Set boundary values
T(1,:) = T_top;                   % Top row
T(end,:) = T_bottom;             % Bottom row
T(:,1) = T_left;                 % Left column
T(:,end) = T_right;              % Right column

% Parameters
eps_s = 1;        % Convergence criterion (1%)
error = 100;      % Initial error
iteration = 0;

% Gauss-Seidel Iteration
while error > eps_s
    error = 0;
    iteration = iteration + 1;
    for i = 2:N+1
        for j = 2:N+1
            prev_iter = T(i,j);
            T(i,j) = 0.25 * (T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1));
            e = abs((T(i,j) - prev_iter) / T(i,j)) * 100;
            if e > error
                error = e;
            end
        end
    end
end

fprintf('Converged in %d iterations with max error = %.4f%%\n', iteration, error);

% Display temperature matrix
disp('Final Temperature Distribution:')
disp(T);

% Plot the result
figure(1)
[X, Y] = meshgrid(0:N+1, 0:N+1);
figure;
surf(X, Y, T);
title('Steady State Temperature Distribution');
xlabel('X'); ylabel('Y'); zlabel('Temperature (°F)');
colorbar;