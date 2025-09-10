% Notes
% Euler's method:
%
%   y_new = y_old + slope * delta_x



% ------------------ P 1.20 ------------------ %
% Problem:  In addition to the downward force of gravity (weight) and drag, an object 
%           falling through a fluid is also subject to a buoyancy force which is 
%           proportional to the displaced volume (Archimedes’ principle). For example, for a sphere
%           with diameter d (m), the sphere’s volume is V = πd^3/6, and its projected area is A = πd^2/4. 
%           The buoyancy force can then be computed as Fb = −ρVg. We neglected buoyancy in our derivation of
%           Eq. (1.8) because it is relatively small for an object like a bungee jumper moving through air. However,
%           for a more dense fluid like water, it becomes more prominent.
%
%           (a) Derive a differential equation in the same fashion as Eq. (1.8), but include the buoyancy force and represent
%           the drag force as described in Sec. 1.4.
%           (b) Rewrite the differential equation from (a) for the special case of a sphere.
%           (c) Use the equation developed in (b) to compute the terminal velocity (i.e., for the steady-state case). 
%           Use the following parameter values for a sphere falling through water: sphere diameter = 1 cm, sphere density = 2700 kg/m3,
%           water density = 1000 kg/m3, and Cd = 0.47.
%           (d) Use Euler’s method with a step size of Δt = 0.03125s to numerically solve for the velocity from t = 0 to 0.25s
%           with an initial velocity of zero.








% ------------------ P 1.23 ------------------ %
% Problem:  As depicted in Fig. P1.23, the downward deflection, y (m), of a cantilever beam 
%           with a uniform load, w = 10,000 kg/m, can be computed as
% 
%           y = (w/24EI)*(x^4 − 4Lx^3 + 6L^2x^2)
% 
%           where x = distance (m), E = the modulus of elasticity = 2 × 10^11 Pa, 
%           I = moment of inertia = 3.25 × 10–4 m4, and L = length = 4 m. This equation can be 
%           differentiated to yield the slope of the downward deflection as a function of x
% 
%           dy/dx = (w/24EI)*(4x^3 − 12Lx^2 + 12L^2x)
% 
%           If y = 0 at x = 0, use this equation with Euler’s method (Δx = 0.125 m) to compute the 
%           deflection from x = 0 to L. Develop a plot of your results along with the analytical
%           solution computed with the first equation.

% Set values
E = 2e11; %Pa
I = 3.25e-4; % m^4
L = 4; % m
w = 10000*9.81; % N/m
delta_x = 0.125;

% x array
x = 0:delta_x:L;
n = length(x);

% y array
y_new = zeros(1,n); % deflection values
y_old = 0;  

% Set slope
slope = @(x) (w/24*E*I)*(4*x^3 - 12*L*x^2 + 12*L^2*x);

% From x = 0 to L
for i = 1:n
    y(i) = y_old + slope(x(i))*delta_x;
    y_old = y(i);
end

plot(x, y_new)





% ------------------ P 4.11 ------------------ %
% Problem:  The Maclaurin series expansion. Starting with the simplest version, cos x = 1, 
%           add terms one at a time to estimate cos(π/3). After each new
%           term is added, compute the true and approximate percent relative errors. 
%           Use your calculator or MATLAB to determine the true value. Add terms until the 
%           absolute value of the approximate error estimate falls below an error criterion 
%           conforming to two significant figures.


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








% ------------------ P 4.12 ------------------ %
% Problem:  Perform the same computation as in Prob. 4.11, but
%           use the Maclaurin series expansion for the sin x to estimate 
%           sin(π/3).


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









% ------------------ P 4.13 ------------------ %
% Problem:  Use zero- through third-order Taylor series expansions
%           to predict f(3) for
%
%           f (x) = 25x^3 − 6x^2 + 7x − 88
%
%           using a base point at x = 1. Compute the true percent relative
%           error for each approximation



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





% ------------------ P 4.16 ------------------ %
% Problem:  Use forward and backward difference approximations of O(h) and a centered difference 
%           approximation of O(h^2) to estimate the first derivative of the function examined in
%           Prob. 4.13. Evaluate the derivative at x = 2 using a step size of h = 0.25. Compare 
%           your results with the true value of the derivative. Interpret your results on the basis
%           of the remainder term of the Taylor series expansion.


% Set values