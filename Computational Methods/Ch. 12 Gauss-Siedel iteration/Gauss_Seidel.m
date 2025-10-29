% ------------------ P 12.1 ------------------ %

% Define matricies
A = [ 3 -0.1 -0.2;
      0.1 7 -0.3;
      0.3 -0.2 10];
B = [7.85;
     -19.3;
     71.4];
x = zeros(1, 3);

% Set C to be have diagonal 0's
C = A;

for i=1:3
    C(i,i) = 0.0;
end

disp(C)

% Implement Gauss Seidel method
max_error = 1;
iter = 0;

while(max_error > 0.10*10^-4)
 
    for i=1:3
        x(i) = (B(i) - C(i, :)*x')/A(i, i);
    end
    error = A * x' - B;
    max_error = max(abs(error));
    iter = iter + 1;

end

disp(x);
disp(iter)



% Can also use Matlab's function fsolve
clear
funcs = @(x) [ 3*x(1) - 0.1*x(2) - 0.2*x(3) - 7.85;
      0.1*x(1) + 7*x(2) - 0.3*x(3) + 19.3;
      0.3*x(1) - 0.2*x(2) + 10*x(3) - 71.4];

X0 = [0; 0; 0];

[x, fvalue] = fsolve(funcs, X0);


% ------------------ P 12.9 ------------------ %

% Additional problem 
% x^2 = 5 - y^2
% y + 1 = x^2

% x^2 -5 + y^2 = 0
% x^2 - 1 - y = 0

clear
funcs = @(x) [ x(1)^2 - 5 + x(2)^2;
               x(2) + 1 - x(1)^2];

% Choose values close to roots
X0 = [0; 0];

[x, fvalue] = fsolve(funcs, X0);
disp(fvalue)