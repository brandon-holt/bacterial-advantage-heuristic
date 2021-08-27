% linear stability analysis
clear; clc; close all;

% system of ODEs and constants
Bmax = 5e5;
Km = 9e12;
a = 1e-13;
b = 2e-2;
r = 3;
kcat = 6e10;
f = @(t,X) [r*X(1)*(1 - (X(1)/Bmax)) - (a*X(1)*X(3));-kcat*X(1)*X(2)/(Km + X(2));(kcat*X(1)*X(2)/(Km + X(2)))- (b*X(1)*X(3))];

% intervals and initial conditions
int = [0 1e4]; %hours
init = [500 4.8e14 0]; 

% range of values: x1 = B, x2 = L, x3 = U
n = 10;
x1 = linspace(0, 1e6, n);
x2 = linspace(0, 1e6, n);
x3 = linspace(0, 1e6, n);
[x,y,z] = meshgrid(x1, x2, x3);

u = zeros(size(x));
v = zeros(size(y));
w = zeros(size(z));
for i = 1:numel(x)
    Xprime = f(int(2), [x(i); y(i); z(i)]);
    u(i) = Xprime(1);
    v(i) = Xprime(2);
    w(i) = Xprime(3);
end

figure;
quiver3(x,y,z,u,v,w); figure(gcf)
xlabel('Bacteria')
ylabel('Locked Drug');
zlabel('Unlocked Drug');