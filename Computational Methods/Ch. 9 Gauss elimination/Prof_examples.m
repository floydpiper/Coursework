% This work is not my own, but rather the work of my professor 

% ------------------ P 9.5 ------------------ %

% We apply the boundary condition at the left and right side
clear;
N = 100; 
T(1)  = 40; T(N) = 200; L = 10; h = 0.01; Ta = 20;
delta_x = L/(N-1); 

% The unknowns are T(2) ... T(N-1). The number of unknowns is N-2

% we create a 1D matrix with N-2 elements for the diagonal terms 
f = (2+h*delta_x^2)*ones(1,N-2);  % diagonal terms
e = -1*ones(1,N-2);               % lower diagonal
g = -1*ones(1,N-2);               % upper diagonal
r = h*delta_x^2*Ta*ones(1,N-2);   % right hand side
r(1) = r(1) + T(1);               % correction for the left B.C.
r(N-2) = r(N-2) + T(N);           % correction for the right B.C.

% Now we call the tridiagonal solver to find the solution
T(2:N-1) = Tridiag(e,f,g,r);

x = 0:delta_x:10;
plot(x,T)
title('Temperature profile'); xlabel('x'); ylabel('Temperature');grid;



% ------------------ P 9.19 ------------------ %

clear;
delta_x = 0.2; E = 250*10^9; I = 3*10^(-4); w = 22500; L = 3;
y(1) = 0; y(16) = 0;
x = 0:delta_x:L;
% we have 14 unknown
f = 2*ones(1,14);
e = -1*ones(1,14);
g = -1*ones(1,14);
r = delta_x^2/(E*I)*(w*x.^2/2 - w*L*x/2);
rr = r(2:15);
y(2:15) = Tridiag(e,f,g,rr); % you solve for your interior points
plot(x,y)
title('Deflection of a beam'); xlabel('x'); ylabel('deflection');grid;

