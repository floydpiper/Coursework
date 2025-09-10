% Set a variable 
a = 2;
b = 5;

% Add variables
c = a + b;

% Multiply variables
d = a * b;

% Subtract variables
e = b - a;

% Divide variables 
f = b / a;

fprintf('%f: ', f)
disp(d);
disp(c);

% Array
x = [1,2,3,4];

% Create array by incrementing - starting at 1 -> 10 by 0.5 
y = 1:0.5:10;

% For each var in y, print i 
for i = y
    fprintf('Current index is %d\n', i);
end

len = length(y);

% Number of elements
fprintf('Number of elements in y is %d\n', len)

% clean output window clc


% Print data 

% %c    Single character
% %d    Decimal notation (signed)
% %e    Exponential notation (using a lowercase e as in 3.1415e+00)
% %E    Exponential notation (using an uppercase E as in 3.1415E+00)
% %f    Fixed-point notation
% %g    The more compact of %e or %f, as defined in [2]. Insignificant zeros do not print.
% %G    Same as %g, but using an uppercase E
% %i    Decimal notation (signed)
% %o    Octal notation (unsigned)
% %s    String of characters
% %u    Decimal notation (unsigned)
% %x    Hexadecimal notation (using lowercase letters a-f)
% %X    Hexadecimal notation (using uppercase letters A-F)

% comment mass CTRL + R

% For loops (more)
for idx = 1:5
    fprintf('Current index is %d\n', idx);
end

for iter = 1:10
    a = a + 1;
    fprintf('New value of a is %d\n', a)
end

% While loops
var = 0;

while var < 10
    disp('Hello')

    var = var + 1;

    % If statement
    if var == 5
        break
    end
end



% Plot values
% @(t) used for functions
g = 9.81; % Acceleration due to gravity
m = 68.1;    % Mass
cd = 0.25; % Drag coefficient


func = @(t)     sqrt(g*m/cd)*tanh(sqrt(g*cd/m)*t)
func(8)

t = 0:0.5:12;

for i = 0:length(t)
    plot(t, func(t));
    xlabel('Time (s)');
    ylabel('Function Value');
    title('Plot of the Function');
end

% Another way 
% 
% for i = 1:length(t)
%     ValSave(i) = func(t(i));
% end
% 
% plot(t, ValSave)
