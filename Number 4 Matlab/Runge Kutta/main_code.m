close all; clear; clc

tend = 10; % Total observation time
ht=0.005; % Time step
n=(tend/ht)+1; % Number of time step

[xSol,ySol] = runkut4(@deqs,0,[0 1],tend,ht);

plot(xSol,ySol)
grid on