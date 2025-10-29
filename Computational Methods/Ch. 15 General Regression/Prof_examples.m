% This work is not my own, but rather the work of my professor 


% ------------------ P 15.1 ------------------ %
clear

x = [0    1    2     3     4     5]';
y = [2.1  7.7  13.6  27.2  40.9  61.1]';
x_sq = x.^2;

n = length(x);

z = [ones(n,1)   x    x_sq]

% The linear equation we are solving is (z'*z)A = z'*y which with 
% the left division can be written as,

A = (z'*z)\(z'*y)

sec_fit = A(1) + A(2)*x + A(3)*x.^2; 

plot(x,y,'og',x,sec_fit);
grid
legend('Data Point','quadratic fit','Location','northwest')
title('Polynomia Regression')
xlabel('X')
ylabel('Y')
 





% ------------------ P 15.6 ------------------ %
clear
U = [0.5 2 10 0.5 2 10 0.5 2 10]';
H = [0.15 0.15 0.15 0.3 0.3 0.3 0.5 0.5 0.5]';
K_L = [0.48 3.9 57 0.85 5 77 0.8 9 92]';
n = length(U);

logK_L = log10(K_L); logU = log10(U); logH = log10(H);

z = [ones(n,1) logU logH];
A = (z'*z)\(z'*logK_L)

a_0 = 10^A(1);

K_L_REG = a_0*U.^A(2).*H.^A(3);

loglog(K_L,K_L_REG,'og',K_L,K_L,'k');
grid
legend('Model Prediction','1:1 line','Location','northwest')
title('Polynomia Regression')
xlabel('K_L Measured')
ylabel('K_L Predicted')
