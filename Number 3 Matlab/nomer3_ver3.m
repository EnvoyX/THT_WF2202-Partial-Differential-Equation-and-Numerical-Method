clc; clear; close all;


% Input SID/NIM
NIM = input("Input NIM/SID: ");
fprintf('SID/NIM Number: %d\n', NIM);
NIM = num2str(NIM);  % Convert to array of string;

% Grid Size (3x3 interior, so 5x5 including boundaries)
N = 5;

% Discretization of Geometry
Length = 5;
width = 5;

% Boundary Conditions from NIM
T_top = str2double(NIM(1:3));
T_right = str2double(NIM(4:5));
T_bottom = str2double(NIM(6:7));
T_left = str2double(NIM(8));

% Convergence criteria
epsilon_s = 1; % convergence tolerance in %
max_iter = 10000;
% Relaxition Paramter for faster calculation
omega = 1; % Relaxation Parameter (1 for Gauss-Seidel)

% Initialize grid
T = zeros(N, N);

% Apply boundary conditions
T(1, :) = T_top;         % Top side
T(N, :) = T_bottom;      % Bottom side
T(:, 1) = T_left;        % Left side
T(:, N) = T_right;       % Right side

% Fix the corners by averaging overlapping BCs
T(1,1)   = (T_top + T_left) / 2;      % Top-left
T(1,N)   = (T_top + T_right) / 2;     % Top-right
T(N,1)   = (T_bottom + T_left) / 2;   % Bottom-left
T(N,N)   = (T_bottom + T_right) / 2;  % Bottom-right

% Initial Matrix T
disp(T);

% Plot initial and boundary conditions
figure;
[x,y] = meshgrid(1:N,1:N);
surf(x,y,T);
colorbar;
xlabel("X (Bottom Plate)"); ylabel("Y (Left Plate)"); zlabel("Temperature (°F)");
title('Initial and Boundary Conditions Temperature Distribution (Gauss-Seidel)');
view(45,30);

% Heatmap IC & BC  Plots
figure;
h = heatmap(T);
h.Title = 'Heated Plate Temperature Gradient Initial & Boundary Conditions';
h.XLabel = 'Left';
h.YLabel = 'Bottom';

% Gauss-Seidel Iteration
for iter = 1:max_iter
    max_error = 0;
    for i = 2:N-1
        for j = 2:N-1
            prev_iter = T(i,j);
            T(i,j) = (1 - omega)*T(i,j) + omega * 0.25 * (T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1));
            error = abs((T(i,j) - prev_iter)/T(i,j)) * 100;
            if error > max_error
                max_error = error;
            end
        end
    end
    
    if max_error < epsilon_s
        fprintf('Converged in %d iterations with max error = %.4f%%\n', iter, max_error);
        break
    end
end

% Display final 3x3 interior values
disp('Final Temperature at Interior Points (3x3):');
disp(T(2:4, 2:4));

% Optional: Full Grid Display
disp('Full Temperature Grid (5x5):');
disp(T);

% Plotting the result
[X, Y] = meshgrid(1:N, 1:N);
figure;
surf(X, Y, T, 'EdgeColor', 'none');
xlabel('X (Bottom Plate)');
ylabel('Y (Left Plate)');
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

