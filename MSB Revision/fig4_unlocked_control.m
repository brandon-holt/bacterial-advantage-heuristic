% the unlocked drug only control for figure 4
close all; clear; clc;

%constants
Bmax = 5e5;
Km = 9e12;
a = 1e-13;
b = 2e-2;

%r and kcat values
R = [0.1, 0.14, 0.3, 0.5, 0.85, 2.2, 1.4, 1.6, 3];
KCAT = [25, 13.6, 22, 15, 13.6, 1, 5.2, 4, 6] * 1e10;

%experimental data
A1 = [336	420	504
100	260	100
28	32	40
0	0	0];

A2 = [420	420	504
160	80	60
16	20	16
120	200	100];

A3 = [336	420	504
20	24	22
2.4	2.4	2.4
0	0	0];

B1 = [336	420	504
20	26	22
2.8	2.2	2.2
0	0	0];

B2 = [336	420	504
26	18	24
3.2	3	2.2
36	18	16];

B3 = [336	420	504
22	34	32
4.2	4	4.4
5.4	4.6	4.4];

C1 = [420	420	504
80	180	160
22	12	16
0	0	0];

C2 = [420	420	504
20	26	24
2.8	3.4	3.4
0	0	0];

C3 = [336	420	504
8	8	12
1.6	3	3
38	30	30];

A1 = [[0; 2; 4; 24], A1];
A2 = [[0; 2; 4; 24], A2];
A3 = [[0; 2; 4; 24], A3];
B1 = [[0; 2; 4; 24], B1];
B2 = [[0; 2; 4; 24], B2];
B3 = [[0; 2; 4; 24], B3];
C1 = [[0; 2; 4; 24], C1];
C2 = [[0; 2; 4; 24], C2];
C3 = [[0; 2; 4; 24], C3];

Data = {A1, A2, A3, B1, B2, B3, C1, C2, C3};
SEE = [];
ModelData = {};

figure();
for i = [1:numel(R)]
    r = R(i);
    kcat = KCAT(i);
    %system of equations
    f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(2)); -(b*x(1)*x(2))];
    %set interval and initial conditions
    int = [0 24]; %hours
    init = [500 4.8e14]; %500 bac/uL (1000X dil), %800 uM drug = 4.8e14 copies/uL
    options1 = odeset('Refine',4);
    options2 = odeset(options1,'NonNegative',1);
    [t,xa] = ode15s(f,int,init,options2);
    
    % plot simulated data
    subplot(3, 3, i);
    plot(t,xa(:,1))
    title('Bacteria')
    xlabel('t'), ylabel('B')
    hold on;
    
    %calculate standard error of the estimate
    %sqrt(sum(Y - Yp)^2/N)
    %sqrt(sum(DIFF^2)/N)
    DIFF = 0;
    N = 0;
    % go through each time point
    [r, c] = size(Data{i});
    for j = [1:r]
        % go through each replicate
        for k = [2:c]
            left = xa(t < Data{i}(j,1), 1);
            right = xa(t > Data{i}(j,1), 1);

            if numel(left) == 0
                left = right(1);
            end
            if numel(right) == 0
                right = left(end);
            end
            
            Y = 0.5 * (left(end) + right(1));
            Yp = Data{i}(j, k);
            
            DIFF = DIFF + ((Y - Yp)^2);
            N = N + 1;
            
            scatter(Data{i}(j,1), Data{i}(j, k));
        end
    end
    
    SEE = [SEE sqrt(DIFF/N)];
    ModelData = [ModelData [t xa(:,1)]];
    
    hold off;
end