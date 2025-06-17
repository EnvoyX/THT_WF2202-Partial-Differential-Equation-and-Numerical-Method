function F = deqs(x,y)
% diiferential equation

F = zeros(1,2);
F(1) = y(2);
F(2)= -5*y(1) - 2.5*y(2);