clear;clc;

% run experiment with default BAH = log10(r/kcat)
i = 1;
timeVersusOutcome_default = {};
figure();
for tTot = [1 2 24 100 500 10000]
    
    outcome_default = [];
    BAH_default = [];
    Bmax = 5e4;
    
    for r = 3*logspace(-2,0,5)
        for kcat = 2.5*logspace(9,11,5)
            for Km = 4.5*logspace(14,15,5)
                for a = logspace(-14, -12, 5)
                    for b = logspace(-3, -1, 5)
                        %system of equations
                        f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(3));-kcat*x(1)*x(2)/(Km + x(2));(kcat*x(1)*x(2)/(Km + x(2)))- (b*x(1)*x(3))];
                        %set interval and initial conditions
                        int = [0 tTot]; %hours, to find steady state value
                        init = [500 2.4e16 0]; %excess substrate
                        options1 = odeset('Refine',4);
                        options2 = odeset(options1,'NonNegative',1);
                        [t,xa] = ode15s(f,int,init, options2);
                        BAH_default = [BAH_default log10((r)/(kcat))];
                        outcome_default = [outcome_default xa(end,1)/Bmax];
                    end
                end
            end
        end
    end
    
    timeVersusOutcome_default = [timeVersusOutcome_default tTot [BAH_default' outcome_default']];
    subplot(2,3,i)
    scatter(BAH_default,outcome_default)
    xlabel('B.A.H.')
    ylabel('bacteria win (1) or lose (0)')
    i = i + 1;
    
end
title('Default BAH = log10(r/kcat)');





% run experiment with modified bah = log10(r*b/kcat*a)
clear;clc;
i = 1;
timeVersusOutcome_modified = {};
figure();
for tTot = [1 2 24 100 500 10000]
    
    outcome_modified = [];
    BAH_modified = [];
    Bmax = 5e4;
    
    for r = 3*logspace(-2,0,5)
        for kcat = 2.5*logspace(9,11,5)
            for Km = 4.5*logspace(14,15,5)
                for a = logspace(-14, -12, 5)
                    for b = logspace(-3, -1, 5)
                        %system of equations
                        f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(3));-kcat*x(1)*x(2)/(Km + x(2));(kcat*x(1)*x(2)/(Km + x(2)))- (b*x(1)*x(3))];
                        %set interval and initial conditions
                        int = [0 tTot]; %hours, to find steady state value
                        init = [500 2.4e16 0]; %excess substrate
                        options1 = odeset('Refine',4);
                        options2 = odeset(options1,'NonNegative',1);
                        [t,xa] = ode15s(f,int,init, options2);
                        BAH_modified = [BAH_modified log10((r * b)/(kcat * a))];

                        outcome_modified = [outcome_modified xa(end,1)/Bmax];
                    end
                end
            end
        end
    end
    
    timeVersusOutcome_modified = [timeVersusOutcome_modified tTot [BAH_modified' outcome_modified']];
    
    subplot(2,3,i)
    scatter(BAH_modified,outcome_modified)
    xlabel('B.A.H.')
    ylabel('bacteria win (1) or lose (0)')
    i = i + 1;
end
title('Modified BAH = log10(r * b / kcat * a)');