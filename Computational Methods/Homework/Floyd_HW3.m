% ------------------ P 6.16 ------------------ %
% Init values
r = 2; 
L = 5; 
V = 8; 

% Define function 
f = @(h) (r^2 * acos((r - h)/r) - (r - h).*sqrt(2*r*h - h.^2)) * L - V;

% Bisection method
h_left = 0;          
h_right = r; 
iter = 0;


while (h_right - h_left)/2 > 1.e-4
    iter = iter + 1;
    h_new = (h_left + h_right)/2;
    if f(h_left) * f(h_new) < 0
        h_right = h_new;
    else
        h_left = h_new;
    end
end


disp("------------------ P 6.16 ------------------")
fprintf('Number of iteration is %d\n', iter)
fprintf('final solution is %3.4f\n', h_new);
fprintf('value of the function is %8.2e\n\n', f(h_new));






% ------------------ P 6.19 ------------------ %
% Init values
R = 225;          
C = 0.6e-6;        
L = 0.5;           
Z_target = 100;    

% Define function
f = @(w) (1 / sqrt((1/R^2) + (w*C - 1./(w*L)).^2)) - Z_target;

w_left = 1;      
w_right = 1000;  


% Bisection method
tol = 1e-6;
iter = 0;

while (abs(w_right - w_left) > tol)
    iter = iter + 1;
    w_new = (w_left + w_right)/2;
    
    if f(w_left)*f(w_new) < 0
        w_right = w_new;
    else
        w_left = w_new;
    end
end

omega = (w_left + w_right)/2;

disp("------------------ P 6.19 ------------------")
fprintf('Angular frequency ω = %.6f rad/s\n', omega);
fprintf('Number of iterations = %d\n', iter);
fprintf('Function value f(ω) = %.6e\n', f(iter));




% ------------------ P 6.24 ------------------ %
% Init values
a = -10;  
b = -6;   
tol = 1e-4;

% Define function
f = @(s) s.^5 + 15*s.^4 + 73*s.^3 + 153*s.^2 + 153*s + 90;

% Bisection method
s_left = a;
s_right = b;
iter = 0;

while (s_right - s_left)/2 > tol
    iter = iter + 1;
    s_new = (s_left + s_right)/2;
    if f(s_left)*f(s_new) < 0
        s_right = s_new;
    else
        s_left = s_new;
    end
end


disp("------------------ P 6.24 ------------------")
fprintf('Number of iterations: %d\n', iter)
fprintf('Final root estimate: %.6f\n', s_new)
fprintf('Function value: %8.2e\n\n', f(s_new))