%% --- MATH135A Project 3: Problem 2 ---

% Compute the root of x - 0.5sin(x) - 1 = 0 with x0_val = 0.5 with
% Secant method with x1 from one iteration of Newton's method
% Fixed point method with choice of g(x)
format long;
% ---- Initial Setup ----
disp(' ');
disp('--- Problem 2 ---');


f2 = @(x) x - 0.5*sin(x) - 1;
df2 = @(x) 1 - 0.5*cos(x); % For Newton's Method again
true_r = 1.498701133517848;
x0 = 0.5;

%% --- PART (a): SECANT METHOD ---
disp('Part (a): Secant Method');

% x1 using one Iteration of Newton's Method
x1 = x0 - f2(x0)/df2(x0);
fprintf('Calculated x1 (via Newton): %.10f\n', x1);

fprintf('%-5s %-15s %-20s %-15s\n', 'Iter', 'x_n', 'f(x_n)', 'Abs Error');

% Initial printing for k=0 and k=1
fprintf('%-5d %-15.10f %-20.10f %-15.10f\n', 0, x0, f2(x0), abs(x0 - true_r));
fprintf('%-5d %-15.10f %-20.10f %-15.10f\n', 1, x1, f2(x1), abs(x1 - true_r));

% Secant Loop (5 iteration)
p_prev = x0;
p_curr = x1;

for k = 2:6
    % Secant Formula
    % x_new = x_curr - f(x_curr) * (x_curr - x_prev) / (f(x_curr) - f(x_prev))

    top = f2(p_curr) * (p_curr - p_prev);
    bot = f2(p_curr) - f2(p_prev);
    p_new = p_curr - (top/bot);

    % Print the current iteration results
    fprintf('%-5d %-15.10f %-20.10f %-15.10f\n', k, p_new, f2(p_new), abs(p_new - true_r));

    % Update previous and current values for the next iteration
    p_prev = p_curr;
    p_curr = p_new;

end

%% --- Part(b): Fixed point method ---

disp(' ');
disp('Part(b): Fixed Point Method');

%We will be using x = g(x), so it will be g(x) = 1 + 0.5sin(x)

g = @(x) 1 + 0.5*sin(x);

x_fp = x0;
fprintf('%-5s %-15s %-15s\n', 'Iter', 'x_n', 'Abs Error');

fprintf('%-5d %-15.10f %-15.10f\n', 0, x_fp, abs(x_fp - true_r));

for k = 1:5
    x_fp = g(x_fp);
    fprintf('%-5d %-15.10f %-15.10f\n', k, x_fp, abs(x_fp - true_r));
end


disp("Does g'(x) converges or diverges?")
dg = @(x) 0.5*cos(x);
if abs(dg(true_r)) < 1
    disp(abs(dg(true_r)));
    disp('Converges');
else
    disp(abs(dg(true_r)));
    disp('Diverges');
end

