% ------------------ P 4.11 % ------------------
% Set values
x = pi/3;
n = 2;
cos_x = 1;
true_val = cos(pi/3);
keepLooping = true;

% Run through first 10 terms of series
while keepLooping

    % Build out Maclaurin series
    term = ((-1)^(n/2)) * (x^n)/factorial(n);
    cos_x = cos_x + term;
    
    % Calculate error
    error = (abs(true_val - cos_x)/abs(true_val)) * 100;

    if error < 0.5
        keepLooping = false;
    end

    % Increment n 
    n = n + 2;

end

% Print results
disp("------------------ P 4.11 ------------------")
fprintf("Estimated value of cos(pi/3): %.4f\n", cos_x)
fprintf("True value of cos(pi/3): %.4f\n", true_val)








% ------------------ P 4.12 % ------------------
% Set values
x = pi/3;
n = 0;
sin_x = 0;
true_val = sin(pi/3);
keepLooping = true;
iter = 0;

% Run through first 10 terms of series
while keepLooping

    % Build out Maclaurin series
    term = ((-1)^n) * (x^(2*n+1)) / factorial(2*n+1);
    sin_x = sin_x + term;
    
    % Calculate error
    error = (abs(true_val - sin_x)/abs(true_val)) * 100;

    if error < 0.5
        keepLooping = false;
    end

    % Increment n 
    n = n + 1;
    
end


% Print results
fprintf("\n\n")
disp("------------------ P 4.12 ------------------")
fprintf("Estimated value of sin(pi/3): %.4f\n", sin_x)
fprintf("True value of sin(pi/3): %.4f\n", true_val)









% ------------------ P 4.13 % ------------------
% Set values
f = @(x) 25*x^3 - 6*x^2 + 7*x - 88;
base_point = 1;
x_val = 3;
true_value = f(x_val);

% Calculate derivatives 
f1 = @(x) 75*x^2 - 12*x + 7;  
f2 = @(x) 150*x - 12;          
f3 = @(x) 150;     

% Taylor series
fprintf("\n\n")
disp("------------------ P 4.13 ------------------")
estimated = f(1);
error = (abs(true_value - estimated)/true_value) * 100;
fprintf("Estimated value of at 0th order: %.4f\n", estimated)
fprintf("Error at 0th order: %.4f\n", error)
fprintf("\n")

estimated = f(1) + f1(1)*2;
error = (abs(true_value - estimated)/true_value) * 100;
fprintf("Estimated value of at 1st order: %.4f\n", estimated)
fprintf("Error at 1st order: %.4f\n", error)
fprintf("\n")

estimated = f(1) + f1(1)*2 + (f2(1)/2)*2^2;
error = (abs(true_value - estimated)/true_value) * 100;
fprintf("Estimated value of at 2nd order: %.4f\n", estimated)
fprintf("Error at 2nd order: %.4f\n", error)
fprintf("\n")

estimated = f(1) + f1(1)*2 + (f2(1)/2)*2^2 + (f3(1)/6)*2^3;
error = (abs(true_value - estimated)/true_value) * 100;
fprintf("Estimated value of at 3rd order: %.4f\n", estimated)
fprintf("Error at 3rd order: %.4f\n", error)





% ------------------ P 4.16 % ------------------
% Set values