% Gauss Seidel Number 1c
clear all; close all; clc
format short


% Input User
n = input('The number of equations/the number of unknown variables: ');
for i=1:n
    for j=1:n
        prompt = sprintf('Enter the element A%d,%d: ', i, j);
        A(i,j)=input(prompt);
    end
end

for i=1:n
    prompt2 = sprintf('Enter the element b%d: ', i);
    b(i,1) = input(prompt2);
end

maxerror = 1e-5;
disp("Matrix A")
disp(A);
disp("Matrix B")
disp(b);

n = length(b);
x = zeros(n,1); % Initial guess
x_prev_iteration = x;


% Gauss-Seidel Formula
for k = 1:7
    for i = 1:n
        sum1 = A(i,1:i-1) * x(1:i-1);
        sum2 = A(i,i+1:n) * x_prev_iteration(i+1:n);
        x(i) = (b(i) - sum1 - sum2)/A(i,i);
    end
x_prev_iteration = x;
end

x_Results = x;

% Display result rounded to 4 decimal places
fprintf('\nGauss-Seidel Iterative Solution after 7 iterations:\n');
for i = 1:n
    val = floor(x_Results(i) * 10000) / 10000; 
    fprintf('X%d = %.4f\n', i, val);
end