% Professors example
% ---------------------------------------------------------------------- %
% Free fall problem
% We want to solve the free fall problem using the Secant method
% Similar to Example 6.7 of the book

clear;
format long
cd = 0.25; g = 9.8; v = 36; t = 4; 

% This defines a function "func" as a function of "m"
func = @(m) sqrt(g*m/cd)*tanh(sqrt(g*cd/m)*t) - v;

delta_x = 1.e-4;
iter = 0;
x_old = 200;
x_new = x_old;

while (abs(func(x_new)) > 1.e-6)

       iter = iter + 1;
       
       x_new =  x_old - ...
             func(x_old)*delta_x/(func(x_old + delta_x) - func(x_old));

       x_old = x_new;    
end   

fprintf('Number of iteration is %d\n', iter)
fprintf('The value of the function is %8.2e\n', func(x_new))
fprintf('The solution is %3.15e\n', x_new)

  % you can also call the Matlab function of fzero

 Sol = fzero(func,200)


% ----------------- Problem 6.20 ----------------- %
% Using the Secant Method

% Init values
k1 = 40000.0;
k2 = 40.0;
m = 95.0;
g = 9.81;
h = 0.43;
range = 0.0:0.01:0.4;
keepLooping = true;

% Define function 
f = @(d) 2*k2*d.^(2.5/5) + 0.5*k1*d.^2 - m*g*d - m*g*h;

% Plot to get value for d to start
for i = 1:length(range)
    Y(i) = f(range(i));
end
plot(range, Y); grid on;

% Set abitrary value for root
d_old = 0.16;
delta_d = 1.e-4;
d_new = d_old;

while (abs(f(d_new))>10^-6)

    d_new = d_old - (f(d_old) * delta_d)/(f(d_old + delta_d) - f(d_old));
    d_old = d_new;
    
end

disp(f(d_new))
disp(d_new)


% Professors example
% ---------------------------------------------------------------------- %


clear;
k1 = 40000; k2 = 40; m = 95; g = 9.81; h = 0.43;
func = @(d) 2*k2*d^2.5/5 + 0.5*k1*d^2 - m*g*d - m*g*h;
dd = 0:0.01:0.4;

for i = 1:length(dd)
    Y(i) = func(dd(i));
end

plot(dd,Y); grid;
delta_x = 1.e-4;
iter = 0;
x_old = 200;
x_new = x_old;

while (abs(func(x_new)) > 1.e-6)
    iter = iter + 1;
    x_new = x_old - ...
    func(x_old)*delta_x/(func(x_old + delta_x) - func(x_old));
    x_old = x_new;
end

fprintf('Number of iteration is %d\n', iter)
fprintf('The value of the function is %8.2e\n', func(x_new))
fprintf('The solution is %3.15e\n', x_new)
% you can also call the Matlab function of fzero
Sol = fzero(func,.1)