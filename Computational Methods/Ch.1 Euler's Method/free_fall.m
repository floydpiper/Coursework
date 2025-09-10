% Using Euler's method to solve the free fall problem

% Init variables 
cd = 0.25; % drift velocity??
g = 9.81; % acceleration due to gravity (m/s^2)
t = 0; % initial time (s)
dt = 0.5; % time step (s)
vNew = 0; % initial velocity (m/s)
vOld  =0;
vals = [0]
t = 0:0.5:12
i = 0;
m = 80; % mass in kg
time = 0.0;

while time < 12
    % Write function
    vNew = vOld + (g - (cd/m) * vOld^2)*dt; 
    % Add velocity to array
    vals = [vals, vNew];
    % Set current value to the old value
    vOld = vNew;
    % Increase the time
    time = time + dt;
end

disp(length(t))
disp(length(vals))
% Plot the function vs time
plot(t, vals, 's--g')
xlabel('Time (s)');
ylabel('Function Value');
title('Plot of the Function');







% Professor's function
% ----------------- %
% Example 1.2, Euler method to solve the free fall problem
clear;
cd = 0.25; g = 9.8; t = 4; m = 80; delta_t = 0.5; t_final = 12;
exact_vel = @(t) sqrt(g*m/cd)*tanh(sqrt(g*cd/m)*t);
time = 0.0; ii = 1;
vel_old = 0.0; % initial condition

while(time<= t_final)
    % store the numerical and exact value in an array to plot them
    vel_num(ii) = vel_old;
    vel_exact(ii) = exact_vel(time);
    % update the velocity using the Euler method
    vel_new = vel_old + (g - cd/m*vel_old^2)*delta_t;
    vel_old = vel_new;
    % update the time
    time = time + delta_t;
    ii = ii + 1;
end

tt = 0.0:delta_t:t_final;
plot(tt,vel_num,'s--g', tt, vel_exact,'-b')
grid
legend('Numerical solution', 'Exact solution','Location','northwest')
title('Plot of velocity in time')
xlabel('Time (second)')
ylabel('Exact and numerical velocities (m/s)')
