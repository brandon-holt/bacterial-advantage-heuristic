clear;clc;
i = 1;
timeVersusOutcome = {};
figure();
for tTot = [1 2 12 24 50 100]
    
    outcome = [];
    BAH = [];
    r_s = [];
    kcat_s = [];
    Bmax_s = [];
    Km_s = [];
    
    for r = 3*logspace(-1,0,10)
        for kcat = 1*logspace(10.5,11.5,10)
            for Bmax = linspace(1e9, 2e9, 5)
                for Km = linspace(1e16, 2e16, 5)
                    %constants
                    a = 3e-16; % mL/hr, +/- 2e-16 (stdev) from 10 min study
                    b = 1e-5; % mL/hr, between 10^-4 and 10^-6 from 10 min study
                    %system of equations
                    f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(3));-kcat*x(1)*x(2)/(Km + x(2));(kcat*x(1)*x(2)/(Km + x(2)))- (b*x(1)*x(3))];
                    %set interval and initial conditions
                    int = [0 tTot]; %hours, to find steady state value
                    init = [4000 7.525e15 0]; %bacteria/mL, 7.525e15 drugs/mL = 12.5 uM drug
                    options1 = odeset('Refine',4);
                    options2 = odeset(options1,'NonNegative',1);
                    [t,xa] = ode15s(f,int,init, options2);
                    BAH = [BAH log10(r/kcat)];
                    r_s = [r_s r];
                    kcat_s = [kcat_s kcat];
                    Bmax_s = [Bmax_s Bmax];
                    Km_s = [Km_s Km];

                    outcome = [outcome xa(end,1)/Bmax];
                end
            end
        end
    end
    
    timeVersusOutcome = [timeVersusOutcome tTot [BAH' outcome']];
    
    subplot(2,3,i)
    scatter(BAH,outcome)
    xlabel('B.A.H.')
    ylabel('bacteria win (1) or lose (0)')
    table = [BAH' outcome' r_s' kcat_s' Bmax_s' Km_s'];
    i = i + 1;
end