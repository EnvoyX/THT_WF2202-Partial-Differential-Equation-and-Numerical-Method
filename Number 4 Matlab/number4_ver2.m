clear; clc; clf;

% Input NIM or SID
NIM_str = input('Enter the SID/NIM number: ');
fprintf('SID/NIM Number: %d\n', NIM_str);

NIM_str = num2str(NIM_str);  % Convert to string for digit access


% Known Data & Constants
global M E I;
M       = str2double(NIM_str(5:8));
E       = str2double(NIM_str(4:8));
I       = 1500;     % Moment of inertia

L = 200;        % mm
x0 = 0;            % Start position or Initial X
y0 = [0, 0];  % y(0) = 0, y'(0) = 0 or Initial Conditions
z0 = 0;
h = 20;   % Step Size

% Solve using Modified Euler (RK2) and RK4 
[x_RK2, y_RK2] = ModifiedEulerRK2(x0, y0, L, h);
[x_RK4, y_RK4] = ClassicalRK4(x0, y0, L, h);

% Plot results
plot(x_RK2, y_RK2(:,1), '-*', 'DisplayName', 'Modified Euler RK2');
hold on;
grid on;
plot(x_RK4, y_RK4(:,1), 'DisplayName', 'Classical Runge-Kutta 4th Order');
xlim([0, L]);
ylim([0, max(y_RK4(:,1))]);
legend('Location', 'best');
xlabel('x (mm)');
ylabel('Deflection y(x) (mm)');
title('Beam Deflection using Bernoulli-Euler Theory');



% Delcare Functions

% ---- ODE System for Beam Deflection ----
function dydx = BeamODE(x, y)
    global M E I;
    slope = y(2);  % dy/dx
    dy1 = slope;
    dy2 = M / (E * I) * (1 + slope^2)^(3/2);  % d²y/dx²
    dydx = [dy1, dy2];
end

% ---- Modified Euler (RK2) Method ----
function [x_values, y_values] = ModifiedEulerRK2(x0, y0, x_end, h)
    numSteps = floor(x_end / h);
    yDim = length(y0);
    x_values = zeros(numSteps+1, 1);
    y_values = zeros(numSteps+1, yDim);
    x_values(1) = x0;
    y_values(1, :) = y0;

    x = x0;
    y = y0;

    for i = 2:numSteps+1
        h = min(h, x_end - x);  % Prevent overshooting
        k1 = BeamODE(x, y);
        k2 = BeamODE(x + h, y + k1 * h);
        y = y + 0.5 * (k1 + k2) * h;
        x = x + h;
        x_values(i) = x;
        y_values(i, :) = y;
    end
end

% ---- Classical Runge-Kutta 4th Order ----
function [x_values, y_values] = ClassicalRK4(x0, y0, x_end, h)
    numSteps = floor(x_end / h);
    yDim = length(y0);
    x_values = zeros(numSteps+1, 1);
    y_values = zeros(numSteps+1, yDim);
    x_values(1) = x0;
    y_values(1, :) = y0;

    x = x0;
    y = y0;

    for i = 2:numSteps+1
        h = min(h, x_end - x);  % Prevent overshooting
        k1 = BeamODE(x, y);
        k2 = BeamODE(x + h/2, y + h/2 * k1);
        k3 = BeamODE(x + h/2, y + h/2 * k2);
        k4 = BeamODE(x + h,   y + h * k3);
        y = y + (k1 + 2*k2 + 2*k3 + k4) * h / 6;
        x = x + h;
        x_values(i) = x;
        y_values(i, :) = y;
    end
end
