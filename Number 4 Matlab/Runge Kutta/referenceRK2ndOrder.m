clear all;
clc;

h2 = 0.5;

% Initial Condition
y3(1) = 4;
dy3 = @(u) u;
u3(1) = 0;
du3 = @(u, y) -(0.6*u) - (8*y);
t = 0:h2:5;
a1 = (t(end)-t(1))/h2;

% Heun Method
for i = 1:a1
    % First Order
    k1y = dy3(u3(i));
    k1u = du3(u3(i), y3(i));
    
    % Second Order
    k2y = dy3(u3(i) + k1u*h2);
    k2u = du3(u3(i) + k1u*h2, y3(i) + k1y*h2);
    
    % Final Condition
    y3(i+1) = y3(i) + ((k1y/2) + (k2y/2))*h2;
    u3(i+1) = u3(i) + ((k1u/2) + (k2u/2))*h2;
end

plot(t, y3), xlabel('t'), ylabel('y_3(t)');
title('Second-order RK Method with h = 0.5')