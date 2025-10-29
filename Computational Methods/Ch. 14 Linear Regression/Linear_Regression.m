% ------------------ P 14.4 ------------------ %
% Create array values for v and F
v = [10 20 30 40 50 60 70 80];
F = [25 70 380 550 610 1220 830 1450];
n = length(v);

% Init values for loop
iter = 0;
sum_x_y = 0;
sum_x = 0;
sum_y = 0;
sum_x_squared = 0;
xy = 0;
x = 0;
y = 0;
x_squared = 0;

% Find the summations beforehand
for i=1:n
    xy = v(i)*F(i);
    sum_x_y = sum_x_y + xy;
end

for i=1:n
    x = v(i);
    sum_x = sum_x + x;
end

for i=1:n
    y = F(i);
    sum_y = sum_y + y;
end

for i=1:n
    x_squared = (v(i))^2;
    sum_x_squared = sum_x_squared + x_squared;
end


% Other way to sum:
%sum_vF = sum(v.*F);   sum_v = sum(v);    sum_F = sum(F);   sum_V2 = sum(v.*v);


% Find a0 and a1 using the equation 
a1 = ((n*sum_x_y) - (sum_x * sum_y))/((n*sum_x_squared) - (sum_x)^2);

average_F = sum_y/n;
average_v = sum_x/n;
a0 = average_F - a1*average_v;

fprintf('a1: %4.2f\n', a1)
fprintf('a0: %4.2f', a0)

regression_line = a0 + a1*v;
plot(v, F, 'og', v, regression_line); grid on;

