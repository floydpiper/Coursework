% ------------------ P 12.11 ------------------ %
% Init values
T_11 = [];
T_21 = [];
T_12 = [];
T_22 = [];

% Write equations for nodes
T_11 = (1/4)*(100 + 75 + T_21 + T_12);
T_21 = (1/4)*(T_11 + 75 + 0 + T_22);
T_12 = (1/4)*(100 + T_11 + T_22 + 25);
T_22 = (1/4)*(T_12 + T_21 + 0 + 25);

% Rearrange equations (move T's to one side)
%   4*T_11 - T_12 - T_21 = 175
%   -T_11 + 4*T_21 - T_22 = 75
%   -T_11 + 4*T_12 - T_22 = 125
%   -T_12 - T_21 + 4*T_22 = 25

% Define matricies 
A = [4 -1 -1  0;
    -1  4  0 -1;
    -1  0  4 -1;
     0 -1 -1  4];

B = [175; 75; 125; 25];

x = zeros(1, 4);

% Set C to be have diagonal 0's
C = A;

for i=1:4
    C(i,i) = 0.0;
end

disp(C)

% Implement Gauss Seidel method
max_error = 1;
iter = 0;

while(max_error > 0.10*10^-4)
 
    for i=1:4
        x(i) = (B(i) - C(i, :)*x')/A(i, i);
    end
    error = A * x' - B;
    max_error = max(abs(error));
    iter = iter + 1;

end

% Check answers 
R = A\B;

% Final variables
T_11 = x(1);
T_21 = x(2);
T_12 = x(3);
T_22 = x(4);

R_11 = R(1);
R_21 = R(2);
R_12 = R(3);
R_22 = R(4);


disp("------------------ P 12.11 ------------------")
disp('Gauss Seidel Method:');
fprintf('T_11 = %.4g\nT_21 = %.4g\nT_12 = %.4g\nT_22 = %.4g\n', x(1), x(2), x(3), x(4));
fprintf('\nVerify results:\n');
fprintf('R_11 = %.4g\nR_21 = %.4g\nR_12 = %.4g\nR_22 = %.4g\n\n', R(1), R(2), R(3), R(4));



% ------------------ P 14.12 ------------------ %
% Create array values for W and A
W = [70 75 77 80 82 84 87 90];
A = [2.10 2.12 2.15 2.20 2.22 2.23 2.26 2.30];
n = length(W);

% Need to linearize to use the power law A = a*W^b
x = log(W);
y = log(A);

% Sum values
sum_xy = sum(x.*y);   
sum_x = sum(x);    
sum_y = sum(y);   
sum_x2 = sum(x.*x);

% Find a0 and a1 using the equation 
a1 = (n*sum_xy - sum_x*sum_y)/(n*sum_x2 - (sum_x)^2);
average_x = sum_x/n;
average_y = sum_y/n;
a0 = average_y - a1*average_x;
regression_line = a0 + a1*x; 

disp("------------------ P 14.12 ------------------")
disp('Linear Regression:');
fprintf('a1: %4.4f\n', a1)
fprintf('a0: %4.4f\n', a0)

% Need to back transform 
plot(W, A, 'og', W, exp(regression_line), '-b'); grid on;
xlabel('W (kg)'); ylabel('A (m^2)');
title('Power-law fit: A = a W^b');
legend('Data','Power-law fit','Location','NorthWest');


% Back-transform parameters and compute prediction
a = exp(a0); b = a1;
Wpred = 95; 
Apred = a * Wpred^b;
fprintf('Predicted A at %d kg: %.3f m^2\n', Wpred, Apred)
fprintf('a: %4.4f\n', a)
fprintf('b: %4.4f\n\n', b)




% ------------------ P 14.29 ------------------ %
% Create array values for W and A
t = [0 5 10 15 20];
p = [100 200 450 950 2000];
n = length(t);

% Need to linearize to use the exponential law A = a*W^b
x = t;
y = log(p);

% Sum values
sum_xy = sum(x.*y);   
sum_x = sum(x);    
sum_y = sum(y);   
sum_x2 = sum(x.*x);

% Find a0 and a1 using the equation 
a1 = (n*sum_xy - sum_x*sum_y)/(n*sum_x2 - (sum_x)^2);
average_x = sum_x/n;
average_y = sum_y/n;
a0 = average_y - a1*average_x;
regression_line = a0 + a1*x; 

disp("------------------ P 14.29 ------------------")
disp('Linear Regression:');
fprintf('a1: %4.4f\n', a1)
fprintf('a0: %4.4f\n', a0)

% Need to back transform 
plot(t, p, 'og', t, exp(regression_line), '-b'); grid on;
xlabel('Years'); ylabel('Number of people');
title('Exponential fit: p = a e^{b t}');
legend('Data','Exponential-law fit','Location','NorthWest');


% Back-transform parameters and compute prediction
a = exp(a0); b = a1;
tpred = 25; 
ppred = a * exp(b*tpred);    
fprintf('Predicted A at in %d years: %.0f people\n', tpred, ppred)
fprintf('a: %4.4f\n', a)
fprintf('b: %4.4f', b)