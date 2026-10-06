function [t, x] = rrk4(f, t0, x0, T, N)
    h = (T - t0) / N; % Step size
    t = linspace(t0, T, N+1); %Time vector
    x = zeros(N+1, 1); % Preallocate solution array
    x(1) = x0;    % Initial condition

    for k = 1:N
        tk = t(k);
        xk = x(k);

        k1 = f(tk, xk);
        k2 = f(tk + h/2, xk + h/2 * k1);
        k3 = f(tk + h/2, xk + h/2 * k2);
        k4 = f(tk + h, xk + h * k3);
        x(k+1) = xk + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end

% --- Problem 1: RK4 Testing ---
f = @(t, x) -2*x + sin(t);
t0 = 0;
T = 10;
x0 = 1;
N_vals = [50, 100, 200, 400];

% Exact reference value at t = 10
x_ref = (1/5)*(2*sin(T) - cos(T)) + (6/5)*exp(-2*T);

fprintf('--- Problem 1 Results ---\n');
fprintf('h\t\t x_N\t\t x_ref\t\t error(h)\t error(h)/error(h/2)\n');
fprintf('----------------------------------------------------------------------\n');

errors = zeros(length(N_vals), 1);
for i = 1:length(N_vals)
    N = N_vals(i);
    h = (T - t0) / N;
    [t, x] = rrk4(f, t0, x0, T, N);
    
    x_N = x(end);
    errors(i) = abs(x_N - x_ref);
    
    if i == 1
        fprintf('%.4f\t %.8f\t %.8f\t %e\t -\n', h, x_N, x_ref, errors(i));
    else
        ratio = errors(i-1) / errors(i);
        fprintf('%.4f\t %.8f\t %.8f\t %e\t %.4f\n', h, x_N, x_ref, errors(i), ratio);
    end
end