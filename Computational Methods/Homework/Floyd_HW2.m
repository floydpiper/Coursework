% ------------------ P 5.15 ------------------ %
% Problem:  Figure P5.15a shows a uniform beam subject to a linearly increasing distributed load. 
%           The equation for the resulting elastic curve is (see Fig. P5.15b) 
% 
%           y =  (w0/(120*E*I*L))*(-x^5 + 2*L^2*x^3 - L^4*x)
% 
%           Use bisection to determine the point of maximum deflection 
%           (i.e., the value of x where dy/dx = 0). Then substitute this value into Eq. (P5.15) 
%           to determine the value of the maximum deflection. Use the following parameter values in your computation: 
% 
%           L = 600 cm, E = 50,000 kN/cm2, I = 30,000 cm4, and 𝑤0 = 2.5 kN/cm.

% Init values
L = 600; % Length in cm
E = 50000; % Elastic modulus in kN/cm^2
I = 30000; % Moment of inertia in cm^4
w0 = 2.5; % Distributed load in kN/cm
dy = @(x) (w0/(120*E*I*L))*(-5*x.^4 + 6*L^2*x.^2 - L^4);

% Plot function fo find init x_left and x_right
x = linspace(0, L, 500);
plot(x, dy(x)); grid on;

% Set x_left and x_right
x_left = 250;
x_right = 300;
x_new = x_right;
iter = 0;

% Start bisection method
while (abs(x_right - x_left) > 1e-4)

    iter = iter + 1;

    % Midpoint
    x_new = (x_left + x_right)/2;   
    
    % Root found
    if dy(x_new) == 0
        break;      
    % Root lies between x_left and x_new
    elseif dy(x_new)*dy(x_left) < 0
        x_right = x_new;
    % Root lies between x_right and x_new
    else
        x_left = x_new;
    end
end


% Print results
disp("------------------ P 5.15 ------------------")
fprintf("Bisection Solution: %.4f (after %d iterations)\n\n", x_new, iter);



% Professors example
% ---------------------------------------------------------------------- %
% Problem 5.15 of Chapra, Bisection method
clear;
L = 600; E = 50000; I = 30000; W0 = 2.5;
% This defines a function "func" as a function of "m"
func = @(x) W0/(120*E*I*L)*(-5*x^4 + 6*L^2*x^2 - L^4);
deflec = @(x) W0/(120*E*I*L)*(-x^5 + 2*L^2*x^3 - L^4*x);
x = 0:1:600;
for i = 1:length(x)
y(i) = func(x(i));
end
plot (x,y); grid;
% Now we start with the iteration of the Secant method.
% First we start with our initial guesses using the "input" command
x_left = 200;
x_right = 300;
% Starting the main while loop
iter = 0;
while (abs(x_right - x_left) > 1.e-4)
iter = iter + 1;
x_new = (x_left + x_right)/2.;
if func(x_new)*func(x_left) > 0 %n_new & x_left are at the same side
x_left = x_new;
else
x_right = x_new;
end
end % end while
fprintf('Number of iteration is %d\n', iter)
fprintf('The final solution is %3.4f, \n' , x_new)
fprintf('The value of the derivative is %8.2e \n', func(x_new));
fprintf('The value of the displacement is %8.2e \n', deflec(x_new));







% ------------------ P 5.20 ------------------ %
% Problem:  For fluid flow in pipes, friction is described by a dimensionless number, 
%           the Fanning friction factor f. The Fanning friction factor is dependent 
%           on a number of parameters related to the size of the pipe and the fluid, 
%           which can all  be  represented by another dimensionless quantity, the
%           Reynolds number Re. A formula that predicts f given Re is the von Karman equation:
%
%           1/sqrt(f) = 4*log10(Re*sqrt(f)) - 0.4;
%
%           Typical values for the Reynolds number for turbulent flow are 10,000 to 500,000 and 
%           for the Fanning friction factor are 0.001 to 0.01. Develop a function that uses bisection to 
%           solve for f given a user-supplied value of Re between 2500 and 1,000,000. Design the 
%           function so that it ensures that the absolute error in the result is Ea,d < 0.000005.

% Init values
Re = 50000; % Reynolds number (user-supplied)
g = @(f) (1./sqrt(f)) - 4*log10(Re*sqrt(f)) + 0.4;

% Plot function to find initial x_left and x_right
f_vals = linspace(0.001, 0.01, 500);
plot(f_vals, g(f_vals)); grid on;


% Set x_left and x_right (bracketing interval)
x_left = 0.001;
x_right = 0.01;
x_new = x_right;
iter = 0;

% Start bisection method
while (abs(x_right - x_left) > 5e-6)

    iter = iter + 1;

    % Midpoint
    x_new = (x_left + x_right)/2;   
    
    % Root found
    if g(x_new) == 0
        break;      
    % Root lies between x_left and x_new
    elseif g(x_new)*g(x_left) < 0
        x_right = x_new;
    % Root lies between x_new and x_right
    else
        x_left = x_new;
    end
end

% Print results
disp("------------------ P 5.20 ------------------")
fprintf("Bisection Solution for f: %.8f (after %d iterations)\n\n", x_new, iter);





% Professors example
% ---------------------------------------------------------------------- %
% Problem 5.20 of Chapra, Bisection method
clear;
Re = 100000;
% This defines a function "func" as a function of "m"
fric = @(f) -1./sqrt(f) + 4*log10(Re*sqrt(f)) - 0.4;
x = 0.0001:0.001:0.02;
for i = 1:length(x)
y(i) = fric(x(i));
end
plot (x,y); grid;
% Now we start with the iteration of the Secant method.
% First we start with our initial guesses using the "input" command
x_left = 0.0001;
x_right = 0.02;
% Starting the main while loop
iter = 0;
while (abs(x_right - x_left) > 5.e-9)
iter = iter + 1;
x_new = (x_left + x_right)/2.;
if fric(x_new)*fric(x_left) > 0 %n_new & x_left are at the same side
x_left = x_new;
else
x_right = x_new;
end
end % end while
fprintf('Number of iteration is %d\n', iter)
fprintf('The value of the function is %8.2e, \n' , fric(x_new));
fprintf('The value of the Darcy friction factor is %8.5e \n', x_new);
