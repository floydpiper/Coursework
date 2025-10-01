% Problem: Estimate the first derivative using truncation and round off
% errors
clear

% set values
fx = @(x) -0.1*x^4 - 0.15*x^3 - 0.5*x^2 - 0.25*x + 1.2;

dfdx = @(x) -0.4*x^3 - 0.45*x^2 - x - 0.25;
x = 0.5;
true_val =  dfdx(0.5);
h = 1;


% Loop through until 10 
for i = 1:10

    % Discretized version of the function
    disc_f = (fx(x + h) - fx(x - h)) / (2*h);
    
    % save the error and h
    error(i) = abs(true_val - disc_f);
    h_save(i) = h;

    % Change the value of h
    h = h / 10;



end

loglog(h_save, error, 's--g')
grid
xlabel('Step size')
ylabel('Error')





% Professors example
% ---------------------------------------------------------------------- %
% Example 4.5 of the book to check the round-off and truncation errors
clear;

dfdx = @(x) -0.4*x^3 - 0.45*x^2 - x - 0.25; 
fx   = @(x) -0.1*x^4 - 0.15*x^3 - 0.5*x^2 - 0.25*x + 1.2;

x_zero = 0.5;
% at x_zero = 0.5 the true value of the derivative is
True_value = dfdx(x_zero);

%The initial value of step-size, delta_h
delta_h = 1.;

for  i = 1:10

     Delta_save(i) = delta_h;
     disc_f =  (fx(x_zero + delta_h) - fx(x_zero - delta_h))/(2*delta_h);
     Error_save(i) = abs( True_value - disc_f);

     delta_h = delta_h/10.;
end

%plot(Delta_save, Error_save,'s--g')
loglog(Delta_save, Error_save,'s--g')
grid
legend('Error','Location','northwest')
title('Plot of error vs. grid size')
xlabel('Grid size')
ylabel('Error')
