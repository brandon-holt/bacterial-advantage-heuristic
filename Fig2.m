clear;clc;
i = 1;
timeVersusOutcome = {};
figure();
for tTot = [1 2 24 100 500 10000]
    
    outcome = [];
    BAH = [];
    r_s = [];
    kcat_s = [];
    Bmax_s = [];
    Km_s = [];
    
    for r = 3*logspace(-2,0,10)
        for kcat = 2.5*logspace(9,11,10)
            for Bmax = 5*logspace(4,5,5)
                for Km = 4.5*logspace(14,15,5)
                    %constants
                    a = 1e-13;
                    b = 2e-2;
                    %system of equations
                    f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(3));-kcat*x(1)*x(2)/(Km + x(2));(kcat*x(1)*x(2)/(Km + x(2)))- (b*x(1)*x(3))];
                    %set interval and initial conditions
                    int = [0 tTot]; %hours, to find steady state value
                    init = [500 2.4e16 0]; %excess substrate
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