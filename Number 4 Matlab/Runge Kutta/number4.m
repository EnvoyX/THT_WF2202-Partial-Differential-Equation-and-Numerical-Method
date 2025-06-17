% MATLAB Code
h = 20;
x = 0:h:20;
y_euler = [0, 0.01772];
z_euler = [0, 0.001772];

y_rk4 = [0, 0.01772];
z_rk4 = [0, 0.001772];

figure;
plot(x, y_euler, 'ro--', x, y_rk4, 'bs--');
legend('Modified Euler y', 'RK4 y');
xlabel('x (mm)');
ylabel('y (deflection)');
title('Comparison of Deflection y(x)');

figure;
plot(x, z_euler, 'ro--', x, z_rk4, 'bs--');
legend('Modified Euler z', 'RK4 z');
xlabel('x (mm)');
ylabel('z = dy/dx');
title('Comparison of Slope z(x)');
