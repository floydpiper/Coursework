% ------------------ P 5.19 ------------------ %
% Problem:  A total charge Q is uniformly distributed around a ring-shaped 
%           conductor with radius a. A charge q is located at a distance x from 
%           the center of the ring (Fig. P5.19). The force exerted on the charge 
%           by the ring is given by
%
%           F - q_small*q_big*x/(4*pi*e_zero*(x^2 + a^2)^1.5)
%
%           where e0 = 8.9 × 10−12 C2/(N m2). Find the distance x where the force 
%           is 1.25 N if q and Q are 2 × 10−5 C for a ring with a radius of 0.85 m

% Problem: Find roots of the following function
% f = @(x) F - (1/(4*pi*e_zero))*((q*Q*x)/((x^2 + a^2)^(1.5)));

% Init values
e_zero = 8.9e-12; % C^2/NM^2
F = 1.25;         % N
q = 2e-5;         % C
Q = q;            % C
a = 0.85;         % radius in m 
x = 0.0:0.1:3;

% Define function 
% Have to add . after x
f = @(x) F - (1/(4*pi*e_zero))*((q*Q*x)./((x.^2 + a.^2).^(1.5)));

% Plot
plot(x, f(x))

% Start bracketing methods 
x_left = 0.2;
x_right = 0.3;
iter = 0;

% Start loop (left root) --- WRONG
while (abs(x_right - x_left)>10^-4)

    iter = iter +1;
    x_new = (x_left + x_right)/2;

    if (f(x_new) + f(x_left) > 0)
        x_left = x_new;
    else
        x_right = x_new;
    end
end

fprintf("First root (leftmost): %.4f\n", x_new)

x_left = 1.2;
x_right = 1.3;
iter = 0;

% Start loop (left root) --- wrong
while (abs(x_right - x_left)>10^-4)

    iter = iter +1;
    x_new = (x_left + x_right)/2;

    if (f(x_new) + f(x_left) > 0)
        x_left = x_new;
    else
        x_right = x_new;
    end
end

fprintf("Second root (rightmost): %.4f\n", x_new)



% Professor's Method
% ---------------------------------------------------------------------- %
x_left = 1;
x_right =2.;
x_new = x_right;
iter = 0;

while(abs(f(x_new))>1.e-4)
    iter = iter +1;
    x_new = x_right - ...
        f(x_right)*(x_left - x_right)/(f(x_left) - f(x_right));

    if (f(x_new) * f(x_left) > 0)
        x_left = x_new;
    else
        x_right = x_new;
    end
end

plot(x,f(x),x_new, f(x_new), '-o'); grid;
fprintf("Solution: %.4f\n", x_new)



% Professors example
% ---------------------------------------------------------------------- %
clear;
e_zero = 8.9*10^(-12); q_small = 2.e-5; q_big = 2.e-5; a = 0.85; F = 1.25;
% This defines a function "force" as a function of "x"
force = @(x) F - q_small*q_big*x/(4*pi*e_zero*(x^2 + a^2)^1.5);
x = 0:0.1:3;

for i = 1:31
    y(i) = force(x(i));
end

plot (x,y); grid;
% The function has one root between 0 and 0.5, and another root between
% 1 and 1.5
x_left = 1;
x_right = 2.;
x_new = x_right;
iter = 0;

while (abs(force(x_new)) > 1.e-4)
    iter = iter + 1;

    % Bisection
    % x_new = (x_left + x_right)/2.;

    % False position
    x_new = x_right - ...
    force(x_right)*(x_left - x_right)/(force(x_left) - force(x_right));
    
    if force(x_new)*force(x_left) > 0
        x_left = x_new;
    else
        x_right = x_new;
    end
end % end while

disp("------------------ P 5.19 ------------------")
fprintf('Number of iteration is %d\n', iter)
fprintf('The final solution is %3.4f, and the value of the function is %8.2e\n', ...
x_new, force(x_new));

