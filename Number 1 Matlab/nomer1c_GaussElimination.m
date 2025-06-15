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

x = zeros(n,1); % Initial guess

% Define function
function [x,det] = gauss(A,b)
% Solves A*x = b by Gauss elimination and computes det(A).
% USAGE: [x,det] = gauss(A,b)

% Gauss Elimination Algorithm
if size(b,2) > 1; b = b'; end
% b must be column vector
n = length(b);


for k = 1:n-1   % Elimination phase
   for i= k+1:n
       if A(i,k) ~= 0
            lambda = A(i,k)/A(k,k);
            A(i,k+1:n) = A(i,k+1:n) - lambda*A(k,k+1:n);
            b(i)= b(i) - lambda*b(k);
        end
    end
end
if nargout == 2; det = prod(diag(A)); end
for k = n:-1:1 % Back substitution phase
    b(k) = (b(k) - A(k,k+1:n)*b(k+1:n))/A(k,k);
end
x = b;
end

x_Results = gauss(A,b);

% Display result rounded to 4 decimal places
fprintf('\nGauss Elimination Iterative Solution after 7 iterations:\n');
for i = 1:n
    val = floor(x_Results(i) * 10000) / 10000; 
    fprintf('X%d = %.4f\n', i, val);
end