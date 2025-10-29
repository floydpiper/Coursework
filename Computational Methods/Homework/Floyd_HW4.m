% ------------------ P 8.10 ------------------ %
% Init values
k1=10; k2=30; k3=30; k4=10;
m1=1;  m2=1;  m3=1;

% Write matricies
A = [ (k1+k2)/m1,   -k2/m1,         0;
      -k2/m2,       (k2+k3)/m2,   -k3/m2;
       0,           -k3/m3,       (k3+k4)/m3 ];

x = [0.05; 
     0.04; 
     0.03]; 

% Solve for acceleration a = -Ax
a = -A*x;
a1 = a(1);
a2 = a(2);
a3 = a(3);

disp("------------------ P 8.10 ------------------")
fprintf('a1: %3.2f m/s^2\n', a1);
fprintf('a2: %3.2f m/s^2\n', a2);
fprintf('a3: %3.2f m/s^2\n\n', a3);




% ------------------ P 9.13 ------------------ %
% Init values
F1=500; F2=1000; 
K=4; N=5; 
yin=0.10; xin=0.00; 

% Set up the matricies
a = F1*ones(N-1,1);                 % sub-diagonal
b = -(F1+F2*K)*ones(N,1);           % main diagonal
c = (F2*K)*ones(N-1,1);             % super-diagonal
A = diag(b) + diag(a,-1) + diag(c,1);

B = zeros(N,1); B(1) = -F1*yin; B(N) = -F2*xin;

% Gaussian elimination 
y = A\B;              
x = K*y;             

disp("------------------ P 9.13 ------------------")
fprintf('y_out = %.8g\nx_out = %.8g\n', y(end), x(1));

