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
function yz__1 = n4ODE(x,y)
    global M E I;
    z = y(1,2);
    y__1 = z;
    z__1 = M/(E*I)*(1+z^2)^(1.5);
    yz__1=[y__1,z__1];
end

% Modified Euler Runge-Kutta 2nd Order
function[xs,ys]=ModifiedEulerRK2(x0,y0,L,h)
    n = floor(L./h);
    order = length(y0);
    xs=zeros(n+1,1);
    ys=zeros(n+1,order);
    xs(1)=x0;
    ys(1,:)=y0;
    i=1;
    x=x0;
    y=y0;
    while x < L
        i=i+1;
        h=min(h,L-x);
        k1 = n4ODE(x,y);
        k2 = n4ODE(x+h,y+k1*h);
        x=x+h;
        y=y+0.5*(k1+k2)*h;
        xs(i)=x;
        ys(i,:)=y;
    end
end

% Classical Runge-Kutta 4th Order
function[xs,ys]=ClassicalRK4(x0,y0,L,h)
    n = floor(L./h);
    order = length(y0);
    xs=zeros(n+1,1);
    ys=zeros(n+1,order);
    xs(1)=x0;
    ys(1,:)=y0;
    i=1;
    x=x0;
    y=y0;
    while x<L
        i=i+1;
        h=min(h,L-x);
        k1 = n4ODE(x,y);
        k2 = n4ODE(x+h/2,y+k1*h/2);
        k3 = n4ODE(x+h/2,y+k2*h/2);
        k4 = n4ODE(x+h,y+k3*h);
        x=x+h;
        y=y+(k1+2*k2+2*k3+k4)*h/6;
        xs(i)=x;
        ys(i,:)=y;
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