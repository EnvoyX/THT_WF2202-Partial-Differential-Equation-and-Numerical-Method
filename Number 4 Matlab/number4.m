clear, clc, clf ;


% SID/NIM input
NIM = input('Enter the SID/NIM number: ');
fprintf('SID/NIM Number: %d\n', NIM);

NIM = num2str(NIM);   % convert to NIM array

% Initialization / Known Information


% Declared E & M by extract it from NIM array
global M E I;
E = str2double(NIM(4:8));
M = str2double(NIM(5:8));

% Contants
I = 1500;
L = 200; % mm

% Boundary Condition
x0 = 0;
y0 = [0,0];
z0 = 0;
h = 20;  % step

% Declare Functions

% ODE System for Beam Deflection
function dydx = beamODE(x, y)
    % Access global constants: Moment, Modulus, Inertia
    global M E I;
    
    % Extract slope (dy/dx) from the second element of y vector
    z = y(1, 2);                      % z = dy/dx
    
    % Compute derivatives
    dy1 = z;                          % dy/dx = z
    dy2 = M / (E * I) * (1 + z^2)^(1.5);  % d²y/dx² = M/(EI)*(1 + (dy/dx)^2)^(3/2)
    
    % Return derivative vector [dy/dx, d²y/dx²]
    dydx = [dy1, dy2];
end

% Modified Euler Runge-Kutta 2nd Order
function [xs, ys] = ModifiedEulerRK2(x0, y0, L, h)
    % Initialize step count based on total length and step size
    n = floor(L / h);                 
    order = length(y0);    % Number of equations in the system (should be 2)
    
    % Initialize arrays for x values and solution y values
    xs = zeros(n+1, 1);            % x values
    ys = zeros(n+1, order);        % y values (each row: [y, dy/dx])
    
    % Set initial values
    xs(1) = x0;                       
    ys(1, :) = y0;
    
    % Initialize variables for iteration
    i = 1;
    x = x0;
    y = y0;
    
    % Begin stepping through the domain until x reaches L
    while x < L
        i = i + 1;                    % Increment index
        h = min(h, L - x);            % Adjust final step size if overshooting
        
        % Calculate slopes at beginning and end of the step
        k1 = beamODE(x, y);             
        k2 = beamODE(x + h, y + k1 * h);
        
        % Update values using Modified Euler method
        y = y + 0.5 * (k1 + k2) * h;
        x = x + h;
        
        % Store current step results
        xs(i) = x;
        ys(i, :) = y;
    end
end

% Classical Runge-Kutta 4th Order
function [xs, ys] = ClassicalRK4(x0, y0, L, h)
    % Initialize step count and number of variables in the system
    n = floor(L / h);                
    order = length(y0);              
    
    % Initialize arrays
    xs = zeros(n+1, 1);              
    ys = zeros(n+1, order);          
    
    % Set initial values
    xs(1) = x0;
    ys(1, :) = y0;
    
    % Initialize variables
    i = 1;
    x = x0;
    y = y0;
    
    % Main iteration loop
    while x < L
        i = i + 1;          % Increment index
        h = min(h, L - x);  % Adjust step size if close to boundary
        
        % Compute intermediate slopes
        k1 = beamODE(x, y);                     
        k2 = beamODE(x + h/2, y + k1 * h/2);
        k3 = beamODE(x + h/2, y + k2 * h/2);
        k4 = beamODE(x + h,   y + k3 * h);
        
        % Update the solution using weighted average of slopes
        y = y + (k1 + 2*k2 + 2*k3 + k4) * h / 6;
        x = x + h;
        
        % Store results
        xs(i) = x;
        ys(i, :) = y;
    end
end




% Call Function
[x2,y2]=ModifiedEulerRK2(x0,y0,L,h);
[x4,y4]=ClassicalRK4(x0,y0,L,h);

plot(x2,y2(:,1),'-*')
hold on
grid on
plot(x4,y4(:,1))
xlim([0 L])
ylim([0 max(y4(:,1))])
legend(["Modified Euler Runge-Kutta 2nd Order"  "Classical Runge-Kutta 4th Order"])
xlabel("x(mm)")
ylabel("y(mm)")
title("Beam Deflection from Bernoulli-Euler Theory")

% Display results in table format
fprintf('\n%-10s %-20s %-20s %-20s %-20s\n', 'x (mm)', 'y_RK2 (mm)', 'dy/dx_RK2', 'y_RK4 (mm)', 'dy/dx_RK4');
fprintf(repmat('-', 1, 90)); fprintf('\n');
for i = 1:length(x2)
    fprintf('%-10.2f %-20.6f %-20.6f %-20.6f %-20.6f\n', ...
        x2(i), y2(i,1), y2(i,2), y4(i,1), y4(i,2));
end