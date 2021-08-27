%constants
Bmax = [50000];
Km = 9e12;
r = 3;
kcat = [6956398505.51782];
a = 1e-6;
b = 2e5;

%system of equations
f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(3));-kcat*x(1)*x(2)/(Km + x(2));(kcat*x(1)*x(2)/(Km + x(2)))- (b*x(1)*x(3))];
%set interval and initial conditions
int = [0 72]; %hours
init = [500 2.4e16 0]; %5000 = 10^9 CFU/mL 1.2e16 = 100 uM drug
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
