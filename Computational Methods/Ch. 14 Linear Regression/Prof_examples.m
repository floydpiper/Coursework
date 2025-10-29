% This work is not my own, but rather the work of my professor 


% ------------------ P 14.6 ------------------ %
clear
v = [10 20 30  40  50  60   70  80];
F = [25 70 380 550 610 1220 830 1450];
n = length(v);

% we can have a linear regression in terms of log(F) and log(v)

v_l = log10(v); F_l = log10(F);

sum_vF = sum(v_l.*F_l); sum_v = sum(v_l);
sum_F = sum(F_l); sum_v2 = sum(v_l.*v_l);

a_1 = (n*sum_vF - sum_v*sum_F)/(n*sum_v2 - sum_v^2)
a_0 = sum_F/n - a_1*sum_v/n;
a_0_l = 10^a_0

Regression_log = a_0_l*v.^a_1;

plot(v,F,'og',v,Regression_log);
grid
legend('Data Point','Linear-log regression','Location','northwest')
title('Linear-log Regression')
xlabel('Velocity')
ylabel('Drag')


