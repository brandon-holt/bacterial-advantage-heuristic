% determine the dependence of other variables (a & b on BAH crit)
close all; clear; clc;

reps = 5;
a_s = zeros(reps^2,1);
b_s = zeros(reps^2,1);
bah_crit_s = zeros(reps^2,1);
i = 1;
for a = 5 * logspace(-15,-11,reps)
    for b = .5 * logspace(-4,0,reps)
        
        tic
        a_s(i) = a;
        b_s(i) = b;
        [bah_crit_s(i), r_limits, kcat_limits] = determine_bah_crit(a,b);
        i = i + 1;
        toc
        
    end
end

scatter(log10(a_s ./ b_s), bah_crit_s);
r_limits_title = sprintf('r: %.1e to %.1e',r_limits(1), r_limits(2));
kcat_limits_title = sprintf('k_{cat}: %.1e to %.1e',kcat_limits(1), kcat_limits(2));
title({'Dependence of BAH_{crit} on a & b'; r_limits_title; kcat_limits_title});
xlabel('log_{10}(a/b)');
ylabel('BAH_{crit}');

function [bah_crit, r_limits, kcat_limits] = determine_bah_crit(a, b)

    r_range = 3*logspace(-1.47,0,10);
    r_limits = [r_range(1), r_range(end)];
    kcat_range = 2.5*logspace(9.45,11,10);
    kcat_limits = [kcat_range(1), kcat_range(end)];
    
    outcome = []; BAH = [];
    for r = r_range
        for kcat = kcat_range
            for Bmax = 5*logspace(4,5,5)
                for Km = 4.5*logspace(14,15,5)
                    %system of equations
                    f = @(t,x) [r*x(1)*(1 - (x(1)/Bmax)) - (a*x(1)*x(3));-kcat*x(1)*x(2)/(Km + x(2));(kcat*x(1)*x(2)/(Km + x(2)))- (b*x(1)*x(3))];
                    %set interval and initial conditions
                    int = [0 1e4]; %hours, to find steady state value
                    init = [500 2.4e16 0]; %excess substrate
                    options1 = odeset('Refine',4);
                    options2 = odeset(options1,'NonNegative',1);
                    [t,xa] = ode15s(f,int,init, options2);
                    BAH = [BAH log10(r/kcat)];
                    outcome = [outcome xa(end,1)/Bmax];
                end
            end
        end
    end
    [BAH, i] = sort(BAH);
    outcome = outcome(i);
    bah_crit = mean([BAH(find(outcome < .1, 1, 'last' )), BAH(find(outcome > .9, 1 ))]);

%     scatter(BAH,outcome)
%     xlabel('B.A.H.')
%     ylabel('bacteria win (1) or lose (0)')
    
end