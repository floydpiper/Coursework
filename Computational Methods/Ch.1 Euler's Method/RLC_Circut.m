% Init variables
iOld = 0;
qOld = 1; 
t = 0:0.01:1; % from t -> 0.1s
R = 200; % Ohms
L = 5; % H
C = 10^-4; % F
iter = 0;
iArray = [0]; % init current
qArray = [1]; % init charge in C
idx = 1;
i = 1;
q = 1;
dt = 0.01;
iter = 1;
time = 0.0; 
t_final = 2;

while (time <= t_final)
    % Save differently 
    i_vals(iter) = iOld;
    q_vals(iter) = qOld;
    tt(iter) = time;

    % Calculate the new current and charge based on the previous values
    iNew = ((-qOld/C*L) - (iOld*R)/L)*dt + iOld;
    qNew = iNew*dt + qOld;

    % Store the current value
    iArray = [iArray, iNew]; 
    qArray = [qArray, qNew]; 

    % Update old current and charge
    iOld = iNew; 
    qOld = qNew;

    time = time + dt;
    iter = iter + 1;
end

disp(length(iArray))
% disp(length(t))
% 
plot(time, i_vals)
plot(time, q_vals)
disp(iArray)
disp(t)




% Professor's function
% ----------------- %
% Problem 1.26, Euler method to solve the RLC

clear;
R = 200; L = 5; C = 0.0001; delta_t = 0.01; t_final = 0.1;

time = 0.0; ii = 1;
i_old = 0.0; q_old = 1;  % initial condition

while(time<= t_final)

    % store the numerical and exact value in an array to plot them    
     i_num(ii) = i_old;
     q_num(ii) = q_old; 
     tt(ii)    = time;

    % update the velocity using the Euler method

     i_new = i_old + (-q_old/(L*C) - i_old*R/L)*delta_t;
     i_old = i_new;

     q_new = q_old + i_old*delta_t;
     q_old = q_new;

    % update the time 
     time = time + delta_t;
     ii = ii + 1;

end

%tt = 0.0:delta_t:t_final;

subplot(1,2,1)
plot(tt,i_num,'s--g')
grid
legend('Current (amp)','Location','northwest')
title('Plot of current in time')
xlabel('Time (second)')
ylabel('Current (amp)')


subplot(1,2,2)
plot(tt, q_num,'-b')
grid
legend('Charge (C)','Location','northeast')
title('Plot of charge in time')
xlabel('Time (second)')
ylabel('Charge (C)')









