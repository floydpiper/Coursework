% Init value
k = 0.19;    % min
Ta = 20;     % C
dt = 2;      % min
Told = 95;   % C
time = 0:dt:20;
Temp = [];
n = length(time);
i = 0;
keepLooping = true;

% Euler's Method
while keepLooping

    Tnew = (-k*(Told - Ta))*dt + Told;
    Temp = [Temp, Tnew]; 
    Told = Tnew;
    i = i + 1;

    if (i == n)
        break;
    end

end

% Plot
plot(time, Temp, 'o-g'); 
hold on;
title('Temperature vs. Time');
xlabel('Time (s)'); 
ylabel('Temperature (C)');
legend('Euler Approximation');
grid on




% -------------------------------------- Part a
% Set values
h = 0.5; %step
x = 0.5;
func = @(x) 5*x^3 - 0.5*x^2 - 6*x + 10;
f_prime = @(x) 15*x^2 - x -6;
true_value = f_prime(x);

% Forward, backward, and centered difference approximations of O(h) and O(h^2)
f_i_prime_forward = (func(x + h) - func(x)) / h;
f_i_prime_backward = (func(x) - func(x - h)) / h;
f_i_prime_centered = (func(x + h) - func(x - h)) / (2*h);

fprintf("\n\n")
disp("------------------ P 4.16 ------------------")
fprintf("True value of f_prime: %.4f\n", true_value)
fprintf("\n")

error = (abs(true_value - f_i_prime_forward)/true_value) * 100;
fprintf("O(h): F_i prime forward difference: %.4f\n", f_i_prime_forward)
fprintf("Error from forward difference: %.4f\n", error)
fprintf("\n")

error = (abs(true_value - f_i_prime_backward)/true_value) * 100;
fprintf("O(h): F_i prime backward difference: %.4f\n", f_i_prime_backward)
fprintf("Error from backward difference: %.4f\n", error)
fprintf("\n")

error = (abs(true_value - f_i_prime_centered)/true_value) * 100;
fprintf("O(h^2): F_i prime centered difference: %.4f\n", f_i_prime_centered)
fprintf("Error from centered difference: %.4f\n", error)
fprintf("\n")