% Problem: Free fall problem to determine mass of bungee jumper using
% bisection

% set values
cd = 0.25; % drift velocity
g = 9.8; % acceleration due to gravity (m/s^2)
t = 4; % initial time (s)
v = 36;
x = 50:1:200;
func = @(m) sqrt((g*m)/cd)*tanh(sqrt((g*cd)/m)*t) - v;


% Loop through eqn
for i = 1:length(x)
    y(i) = func(x(i));
end

% Plot the graph to determine what values of x_left and x_right to use
plot(x,y)
grid

% Find the root
x_left = 100;
x_right = 200;
iter = 0;

% Find value of the mass: mass is when x = 0
while(abs(x_right - x_left) > 10e-4)
    iter = iter + 1;

    % This is the bisection method 
    x_new = (x_left + x_right) / 2.;

    if func(x_new) * func(x_left) > 0
        x_left = x_new;
    else
        x_right = x_new;
    end

end

disp(x_new)






% Professors example
% ---------------------------------------------------------------------- %
clear;

% set your parameters

cd = 0.25; g = 9.8; v = 36; t = 4; 

% define your function

func = @(m)  sqrt(g*m/cd)*tanh(sqrt(g*cd/m)*t) - v; 

% If you like you can plot the function first

% I define an array of x like this

x = 50:1:200; 
for i = 1:length(x)
    y(i)  = func(x(i));
end     

plot(x,y); grid;

% choose two values for x_u and x_l

x_left = 100;
x_right = 200;

%  Have a while loop, and check for the sign of your function evaluated at
%  xr = (xu + xl)/2

iter = 0;
while (abs(x_right - x_left)  > 1.e-4)
   
    iter = iter + 1;
    x_new = (x_left+x_right)/2.;

      if  func(x_new)*func(x_left) > 0 
         x_left = x_new;
      else 
         x_right = x_new;
      end 
end 

% you can use fprintf to display some text along with your values
% you will use %d for integer,
% you can use %3.4f for real numbers
% you can use %8.2e for smaller numbers

fprintf('Number of iteration is %d\n', iter)
fprintf('final solution is %3.4f\n,', x_new);
fprintf('value of the function is %8.2e\n,', func(x_new));

% What we learned
% 1. how to use if statement
% 2. how to fprintf for integer, real values and smaller values












