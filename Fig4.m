%constants
Bmax = 5e5;
Km = 9e12;
a = 1e-13;
b = 2e-2;

%r and kcat values
R = [0.1, 0.14, 0.3, 0.5, 0.85, 2.2, 1.4, 1.6, 3];
KCAT = [25, 13.6, 22, 15, 13.6, 1, 5.2, 4, 6] * 1e10;

%experimental data
A1 = [0.000000	500.	500.	500.
1.500000	500.	500.	600.
3.000000	1000.	0.	0.
24.000000	0.	0.	0.];

A2 = [0.000000	500.	500.	500.
4.000000	470.	230.	180.
24.000000	0.	0.	0.];

A3 = [0.000000	500.	500.	500.
2.000000	100.	60.	60.
4.000000	10.	10.	40.
24.000000	0.	0.	0.];

B1 = [0.000000	500.	500.	500.
4.000000	160.	300.	300.
8.000000	230.	90.	400.
24.000000	0.	0.	0.];

B2 = [0.000000	500.	500.	500.
4.000000	1400.	5000.	3000.
8.000000	12000.	9000.	6000.
24.000000	500000.	500000.	500000.];

B3 = [0.000000	500.	500.	500.
2.000000	21000.	24000.	22500
4.000000	500000.	500000.	500000.
24.000000	500000.	500000.	500000.];

C1 = [0.000000	500.	500.	500.
2.500000	500000.	500000.	500000.
5.000000	500000.	500000.	500000.
24.000000	500000.	500000.	500000.];

C2 = [0.000000	500.	500.	500.
2.000000	20000.	500000.	500000.
4.000000	500000.	500000.	500000.
24.000000	500000.	500000.	500000.];

C3 = [0.000000	500.	500.	500.
1.500000	20000.	22000.	31000.
3.000000	500000.	500000.	500000.
24.000000	500000.	500000.	500000.];


Data = {A1, A2, A3, B1, B2, B3, C1, C2, C3};
SEE = [];
ModelData = {};

figure();
for i = [1:numel(R)]
    r = R(i);
    kcat = KCAT(i);
    %system of equations
    f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(3));-kcat*x(1)*x(2)/(Km + x(2));(kcat*x(1)*x(2)/(Km + x(2)))- (b*x(1)*x(3))];
    %set interval and initial conditions
    int = [0 24]; %hours
    init = [500 4.8e14 0]; %500 bac/uL (1000X dil), %800 uM drug = 4.8e14 copies/uL
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