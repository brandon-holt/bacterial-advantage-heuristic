%% experimental data
for i = 1
Time = [0.00 3.25 6.50 22.25]';
NoDrug0 = [3400.000000	4200.000000	4600.000000
3400.000000	3600.000000	3000.000000
3800.000000	4800.000000	4000.000000
36000.000000	54000.000000	54000.000000];
NoDrug75 = [3400.000000	4200.000000	4600.000000
520000.000000	940000.000000	640000.000000
2.800000e+008	5.600000e+008	2.800000e+008
1.980000e+009	1.200000e+009	1.180000e+009];
FreeDrug0 = [3400.000000	4200.000000	4600.000000
0.000000	0.000000	0.000000
1000.000000	800.000000	1200.000000
0.000000	0.000000	0.000000];
FreeDrug75 = [3400.000000	4200.000000	4600.000000
0.000000	0.000000	0.000000
0.000000	0.000000	0.000000
0.000000	0.000000	0.000000];
ProDrug0 = [3400.000000	4200.000000	4600.000000
2600.000000	3200.000000	2400.000000
3200.000000	2400.000000	3400.000000
200.000000	0.000000	200.000000];
ProDrug75 = [3400.000000	4200.000000	4600.000000
42000.000000	36000.000000	40000.000000
64000.000000	68000.000000	56000.000000
28000.000000	28000.000000	32000.000000];
end
figure();
%% Parameters constant

modelData = {'No Drug, 75%', 'No Drug, 0%', 'Free Drug, 75%', 'Free Drug, 0%', 'Pro Drug, 75%', 'Pro Drug, 0%'};

Bmax = 1.5e9;
a = 3e-16; % mL/hr, +/- 2e-16 (stdev) from 10 min study
b = .33e-6; % mL/hr, between 10^-4 and 10^-6 from 10 min study
% b = 0.35e-6; % if 99.9% of drug is left after 10 mins
kcat = 1e11; % h^-1, 0.1 mM TM-TMP by 10 mM GSH in 10 mins, estimate concentration of thiol-reducing agents per unit volume of bacteria 1e10 to 1e12
Km = 1e16; % mL^-1, from 20 uM uptake in 2 hours (Murthy)


%% No Drug Model
% No drug model, 75 %
%system of equations
r = 1.7;
f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax))];
%set interval and initial conditions
int = [0 25]; %hours
init = [4000]; %bacteria/mL, 1.2e16 = 100 uM drug
options1 = odeset('Refine',4);
options2 = odeset(options1,'NonNegative',1);
[t,xa] = ode15s(f,int,init,options2);
%plot results
subplot(2,3,1)
plot(t,xa(:,1))
set(gca, 'YScale', 'log')
hold on;
for i = 1:3
    scatter(Time, NoDrug75(:,i));
    hold on;
end
title('No Drug, 75% Broth')
xlabel('t'), ylabel('B')

modelData{2, 1} = [t,xa(:,1)];

% No drug model, 0 %
%system of equations
r = 0.1;
f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax))];
%set interval and initial conditions
int = [0 25]; %hours
init = [4000]; %bacteria/mL, 1.2e16 = 100 uM drug
options1 = odeset('Refine',4);
options2 = odeset(options1,'NonNegative',1);
[t,xa] = ode15s(f,int,init,options2);
%plot results
subplot(2,3,4)
plot(t,xa(:,1))
hold on;
for i = 1:3
    scatter(Time, NoDrug0(:,i));
    hold on;
end
title('No Drug, 0% Broth')
xlabel('t'), ylabel('B')

modelData{2, 2} = [t,xa(:,1)];

%% Free Drug Model
%Free drug model, 75 %
%system of equations
r = 1.7;
f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(2)); -(b*x(1)*x(2))];
%set interval and initial conditions
int = [0 25]; %hours
init = [4000 7.525e15]; %bacteria/mL, 7.525e15 drugs/mL = 12.5 uM drug
options1 = odeset('Refine',4);
options2 = odeset(options1,'NonNegative',1);
[t,xa] = ode15s(f,int,init,options2);
%plot results
subplot(2,3,2)
plot(t,xa(:,1))
hold on;
for i = 1:3
    scatter(Time, FreeDrug75(:,i));
    hold on;
end
title('Free Drug, 75% Broth')
xlabel('t'), ylabel('B')

modelData{2, 3} = [t,xa(:,1)];

%Free drug model, 0 %
%system of equations
r = 0.1;
f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(2)); -(b*x(1)*x(2))];
%set interval and initial conditions
int = [0 25]; %hours
init = [4000 7.525e15]; %bacteria/mL, 7.525e15 drugs/mL = 12.5 uM drug
options1 = odeset('Refine',4);
options2 = odeset(options1,'NonNegative',1);
[t,xa] = ode15s(f,int,init,options2);
%plot results
subplot(2,3,5)
plot(t,xa(:,1))
hold on;
for i = 1:3
    scatter(Time, FreeDrug0(:,i));
    hold on;
end
title('Free Drug, 0% Broth')
xlabel('t'), ylabel('B')

modelData{2, 4} = [t,xa(:,1)];

%% Pro Drug Model
%Pro drug model, 75 %
%system of equations
r = 1.7;
f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(3));-kcat*x(1)*x(2)/(Km + x(2));(kcat*x(1)*x(2)/(Km + x(2)))- (b*x(1)*x(3))];
%set interval and initial conditions
int = [0 25]; %hours
init = [4000 7.525e15 0]; %bacteria/mL, 7.525e15 drugs/mL = 12.5 uM drug
options1 = odeset('Refine',4);
options2 = odeset(options1,'NonNegative',1);
[t,xa] = ode15s(f,int,init,options2);
%plot results
subplot(2,3,3)
plot(t,xa(:,1))
set(gca, 'YScale', 'log')
hold on;
for i = 1:3
    scatter(Time, ProDrug75(:,i));
    hold on;
end
title('Pro Drug, 75% Broth')
xlabel('t'), ylabel('B')

modelData{2, 5} = [t,xa(:,1)];


%Pro drug model, 0 %
%system of equations
r = 0.1;
f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(3));-kcat*x(1)*x(2)/(Km + x(2));(kcat*x(1)*x(2)/(Km + x(2)))- (b*x(1)*x(3))];
%set interval and initial conditions
int = [0 25]; %hours
init = [4000 7.525e15 0]; %bacteria/mL, 7.525e15 drugs/mL = 12.5 uM drug
options1 = odeset('Refine',4);
options2 = odeset(options1,'NonNegative',1);
[t,xa] = ode15s(f,int,init,options2);
%plot results
subplot(2,3,6)
plot(t,xa(:,1))
hold on;
for i = 1:3
    scatter(Time, ProDrug0(:,i));
    hold on;
end
title('Pro Drug, 0% Broth')
xlabel('t'), ylabel('B')

modelData{2, 6} = [t,xa(:,1)];

