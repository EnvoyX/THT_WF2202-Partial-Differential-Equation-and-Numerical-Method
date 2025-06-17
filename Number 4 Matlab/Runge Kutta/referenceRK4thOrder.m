clear all;
clc;

h3 = 0.5;

% Initial Condition
y4(1) = 4;
dy = @(u) u;
u(1) = 0;
du = @(u,y) -0.6*u - 8*y;
t = 0:h3:5;
a = (t(end)-t(1))/h3;

% Fourth-Order RK Calculation
for i = 1:a
    % First Order
    k1_y(i) = u(i);
    k1_u(i) = -0.6*k1_y(i) - 0.8*y4(i);
    
    % Second Order
    k2_y(i) = u(i) + h3*k1_u(i)/2;
    k2_u(i) = -0.6*k2_y(i) - 0.8*(y4(i) + h3*k1_y(i)/2);
    
    % Third Order
    k3_y(i) = u(i) + h3*k2_u(i)/2;
    k3_u(i) = -0.6*k3_y(i) - 0.8*(y4(i) + h3*k2_y(i)/2);
    
    % Fourth Order
    k4_y(i) = u(i) + h3*k3_u(i);
    k4_u(i) = -0.6*k4_y(i) - 0.8*(y4(i) + h3*k3_y(i));
    
    % Final Condition
    y4(i+1) = y4(i) + 0.5 * (k1_y(i) + 2*k2_y(i) + 2*k3_y(i) + k4_y(i)) / 6;
    u(i+1) = u(i) + 0.5 * (k1_u(i) + 2*k2_u(i) + 2*k3_u(i) + k4_u(i)) / 6;
end

plot(t, y4), xlabel('t'), ylabel('y_4(t)');
title('Fourth-order RK method with h = 0.5')