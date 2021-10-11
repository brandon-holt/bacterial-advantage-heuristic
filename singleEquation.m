%constants
Bmax = 5e5;
Km = 9e12;
a = 1e-13;
b = 2e-2;
r = 3;
kcat = 25e10;

%system of equations
f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(3));-kcat*x(1)*x(2)/(Km + x(2));(kcat*x(1)*x(2)/(Km + x(2)))- (b*x(1)*x(3))];
%set interval and initial conditions
int = [0 24]; %hours
init = [500 4.8e14 0]; %500 bac/uL (1000X dil), %800 uM drug = 4.8e14 copies/uL
options1 = odeset('Refine',4);
options2 = odeset(options1,'NonNegative',1);
[t,xa] = ode15s(f,int,init,options2);

%plot results
figure
subplot(3,1,1)
plot(t,xa(:,1))
title('Bacteria')
xlabel('t'), ylabel('B')
hold on
subplot(3,1,2)
plot(t,xa(:,2))
title('Locked Drug')
xlabel('t'), ylabel(('LD'));
hold on
subplot(3,1,3)
plot(t,xa(:,3))
title('Unlocked Drug')
xlabel('t'), ylabel('UD')
