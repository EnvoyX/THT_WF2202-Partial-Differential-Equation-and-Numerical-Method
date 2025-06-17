function [xSol,ySol] = runkut4(deqs,x,y,xstop,h)
% 4th order runge kutta
% usage: [xXol,ySol] = runkut4(deqs,x,y,xstop,h)
% input:
% deqs = 1st order diff eq
%        F(x,y) = [dy1/dx dy2/dx ...]
% x,y = initial values; y must be row vector y = [a b c d]
% xstop = terminal value of x
% h = increment of x used in integration
% output:
% xsol = x values which solution is computed
% ysol = values of y corresponding to the x-values

if size(y,1) > 1
    y = y'; % y must be row vector
end
xSol = zeros(2,1);
ySol = zeros(2,length(y));
xSol(1) = x;
ySol(1,:) = y;

i = 1;
while x < xstop
    i = i + 1;
    h = min(h,xstop - x);
    k1 = h*feval(deqs,x,y);
    k2 = h*feval(deqs,x + h/2,y + k1/2);
    k3 = h*feval(deqs,x + h/2,y + k2/2);
    k4 = h*feval(deqs,x + h,y + k3);
    y = y + (k1 + 2*k2 + 2*k3 + k4)/6;
    x = x + h;
    xSol(i) = x;
    ySol(i,:) = y; % store current sol
end