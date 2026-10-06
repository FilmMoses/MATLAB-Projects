%% --- MATH135A Project 3: Problem 1 ---
% Compute the zeros of f(x) accuracy to at least 10^(-6) with Bisection and
% newton methods

format long;
f = @(x) exp(sin(x)) + x.^4 - 2*x.^3 - x.^2 -1; %original
df = @(x) cos(x).*exp(sin(x)) + 4*x.^3 - 6*x.^2 - 2*x; %derivative for Newton

true_root_neg =  -0.714889460738844;
true_root_pos = 0.678512694194867;
true_root_large = 2.345483689913113;

%% ----- Part(a): Bisection method ------
disp('----- Problem 1(a): Bisection Method ---');

tol = 1e-6;

a1 = -2; b1 = -0.5;

a2 = 0.5; b2 = 2;

n_step1 = ceil(log2((b1-a1)/tol));
n_step2 = ceil(log2((b2-a2)/tol));

fprintf('Theoretical steps needed for interval 1: %d\n', n_step1);
fprintf('Theoretical steps needed for interval 2: %d\n\n', n_step2);

run_bisection(f, a1, b1, true_root_neg, 'Interval [-2, -0.5]');
run_bisection(f, a2, b1, true_root_pos, 'Interval [0.5, 2]');


%% ----- Part(b): Newton's method ------
disp("--- Problem 1(b): Newtons's Method");

x0_vals = [-0.43, 3];

for k = 1:length(x0_vals)
    x_curr = x0_vals(k);
    fprintf('\nNewton Method with x0 = %.2f\n', x_curr);
    fprintf('%-5s %-15s %-20s %-15s\n', 'Iter', 'x_n', 'f(x_n)', 'Abs Error');
    iter = 0;
    err = 1; % initial error

    while err > 1e-6
        
        %Newton step
        x_next = x_curr - f(x_curr)/df(x_curr);
        %calc error
        err = abs(x_next - true_root_large);

        %update
        % Update current value
        x_curr = x_next;
        iter = iter + 1;
        
        % Print current iteration results
        fprintf('%-5d %-15.8f %-20.8f %-15.8f\n', iter, x_curr, f(x_curr), err);

        if iter > 20, break; end
    end
end


% --- HELPER FUNCTION FOR BISECTION ---
function run_bisection(f, a, b, true_root, label)
    fprintf('Bisection Results for %s:\n', label);
    fprintf('%-5s %-15s %-20s %-15s\n', 'Iter', 'c', 'f(c)', 'Abs Error');

    for k = 1:5
        c = (a + b) / 2;
        fc = f(c);
        err = abs(c - true_root);
        fprintf('%-5d %-15.10f %-20.10f %-15.10f\n', k, c, fc, err);

        if f(a) * fc < 0
            b = c;
        else
            a = c;
        end
    end
    fprintf('\n');
end


