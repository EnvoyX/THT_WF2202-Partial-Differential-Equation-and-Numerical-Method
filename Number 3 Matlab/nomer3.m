clear all; close all ; clc


% Input SID/NIM
NIM = input("Input NIM/SID: ");
fprintf('SID/NIM Number: %d\n', NIM);

NIM = num2str(NIM);  % Convert to array of string;

% Discretization of Geometry
Length = 5;
Width = 5;

elX = 10;% number of elements in x direction
elY = 10; % number of elements in y direction
Nx = elX + 1; % number of nodes in x direction
Ny = elY + 1; % number of nodes in y direction

dx = linspace(0, Length, Nx);
dy = linspace(0, Width, Ny);

% Boundary Condition & Initial Condition
T = zeros(Nx, Ny); % Initialize solution for T

% Boundary Conditions from NIM
T_top = str2double(NIM(1:3));
T_right = str2double(NIM(4:5));
T_bottom = str2double(NIM(6:7));
T_left = str2double(NIM(8));
fprintf('T_top = %d, T_right = %d, T_bottom = %d, T_left = %d\n', T_top, T_right, T_bottom, T_left);


% Apply boundary conditions
T(1, :) = T_top;         % Top side
T(Nx, :) = T_bottom;      % Bottom side
T(:, 1) = T_left;        % Left side
T(:, Ny) = T_right;       % Right side

% Fix the corners by averaging overlapping BCs
T(1,1)   = (T_top + T_left) / 2;      % Top-left
T(1,Ny)   = (T_top + T_right) / 2;     % Top-right
T(Nx,1)   = (T_bottom + T_left) / 2;   % Bottom-left
T(Nx,Ny)   = (T_bottom + T_right) / 2;  % Bottom-right

fprintf(' Top_side = %d \n Bottom_side = %d \n Left_side = %d \n Right_side = %d\n', T(1,1), T(1,Ny), T(Nx,1), T(Nx,Ny));

% Initial Matrix T
disp(T);

% PLot initial and boundary conditions
figure;
[x,y] = meshgrid(dx,dy);
surf(x,y,T);
colorbar;
xlabel("X (Bottom Plate)"); ylabel("Y (Right Plate)"); zlabel("Temperature (°F)");
title('Initial and Boundary Conditions Temperature Distribution (Gauss-Seidel)');
view(45,30);

% Heatmap IC & BC  Plots
figure;
h = heatmap(T);
h.Title = 'Heated Plate Temperature Gradient Initial & Boundary Conditions';
h.XLabel = 'X (Bottom)';
h.YLabel = 'Y (Left)';


% Numerical procedure
% Convergence criteria
epsilon_s = 1; % convergence tolerance in %
max_iter = 10000; % maximum iteration 

% Relaxition Paramter for faster calculation
omega = 1; % Relaxation Parameter (1 for Gauss-Seidel)

% Gauss-Seidel Algorithm 
for iter = 1:max_iter
    max_error = 0;
    for i = 2:Nx-1
        for j = 2:Ny-1
            T_prev_iter = T(i,j);
            T(i,j) = (1 - omega)*T(i,j) + omega * 0.25 * (T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1));
            error = abs((T(i,j) - T_prev_iter)/T(i,j)) * 100;
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

% Display final interior values
fprintf('Final Temperature at Interior Points (%.f x %.f):\n',elX-1,elY-1);
disp(T(2:Nx-1, 2:Ny-1));

% Full Grid Display
fprintf('Full Temperature Grid (%.f x %.f):\n',elX,elY);
disp(T);

% Plotting the result
[X, Y] = meshgrid(dx, dy);
figure;
surf(X, Y, T, 'EdgeColor', 'none');
xlabel('X (Bottom Plate)');
ylabel('Y (Right Plate)');
zlabel('Temperature (°F)');
title('Steady-State Temperature Distribution (Gauss-Seidel)');
colorbar;
view(45,30);


% Heatmap Plots
figure;
h = heatmap(T);
h.Title = 'Heated Plate Temperature Gradient';
h.XLabel = 'X (Bottom)';
h.YLabel = 'Y (Left)';

