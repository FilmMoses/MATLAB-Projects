% --- Problem 2: Linear Solvers ---
n_vals = [10, 100];
beta = 0.1;
max_iter = 10000;
tol = 1e-10;

for idx = 1:length(n_vals)
    n = n_vals(idx);
    
    % Construct Matrix A and vector b
    d = (1:n)';
    A = diag(d) + beta * ones(n, n);
    
    % Ground truth solution
    x_star = zeros(n, 1);
    for i = 1:n
        x_star(i) = sin(i*pi / (n+1)) + i/n;
    end
    b = A * x_star;
    
    fprintf('\n--- Results for n = %d ---\n', n);
    
    % 1. Cholesky
    tic;
    x_cho = my_cholesky_solve(A, b);
    t_cho = toc;
    res_cho = norm(A*x_cho - b) / norm(b);
    fprintf('Cholesky:\t Time = %.6f s,\t Rel. Residual = %e\n', t_cho, res_cho);
    
    % 2. Gauss-Seidel
    tic;
    [x_gs, iter_gs] = my_gauss_seidel(A, b, tol, max_iter);
    t_gs = toc;
    res_gs = norm(A*x_gs - b) / norm(b);
    diff_gs = norm(x_gs - x_cho);
    fprintf('Gauss-Seidel:\t Time = %.6f s,\t Iterations = %d,\t Rel. Res = %e,\t Diff from Cho = %e\n', ...
            t_gs, iter_gs, res_gs, diff_gs);
            
    % 3. Conjugate Gradient
    tic;
    [x_cg, iter_cg] = my_conjugate_gradient(A, b, tol, max_iter);
    t_cg = toc;
    res_cg = norm(A*x_cg - b) / norm(b);
    diff_cg = norm(x_cg - x_cho);
    fprintf('CG:\t\t Time = %.6f s,\t Iterations = %d,\t Rel. Res = %e,\t Diff from Cho = %e\n', ...
            t_cg, iter_cg, res_cg, diff_cg);
end

% Functions

function x = my_cholesky_solve(A, b)
    n = size(A, 1);
    L = zeros(n, n);
    % Cholesky Factorization A = L*L^T
    for j = 1:n
        for i = j:n
            if i == j
                L(i, j) = sqrt(A(i, j) - sum(L(i, 1:j-1).^2));
            else
                L(i, j) = (A(i, j) - sum(L(i, 1:j-1) .* L(j, 1:j-1))) / L(j, j);
            end
        end
    end
    % Forward substitution L*y = b
    y = zeros(n, 1);
    for i = 1:n
        y(i) = (b(i) - L(i, 1:i-1)*y(1:i-1)) / L(i, i);
    end
    % Backward substitution L^T*x = y
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - L(i+1:n, i)'*x(i+1:n)) / L(i, i);
    end
end

function [x, iter] = my_gauss_seidel(A, b, tol, max_iter)
    n = size(A, 1);
    x = zeros(n, 1); % Initial guess x_0 = 0
    norm_b = norm(b);
    
    for iter = 1:max_iter
        x_old = x;
        for i = 1:n
            sum1 = A(i, 1:i-1) * x(1:i-1);
            sum2 = A(i, i+1:n) * x_old(i+1:n);
            x(i) = (b(i) - sum1 - sum2) / A(i, i);
        end
        
        r = b - A*x;
        if (norm(r) / norm_b) < tol
            break;
        end
    end
end

function [x, iter] = my_conjugate_gradient(A, b, tol, max_iter)
    n = size(A, 1);
    x = zeros(n, 1); % Initial guess x_0 = 0
    r = b - A*x;
    p = r;
    norm_b = norm(b);
    
    for iter = 1:max_iter
        Ap = A*p;
        alpha = (r'*r) / (p'*Ap);
        x = x + alpha*p;
        r_new = r - alpha*Ap;
        
        if (norm(r_new) / norm_b) < tol
            break;
        end
        
        beta_val = (r_new'*r_new) / (r'*r);
        p = r_new + beta_val*p;
        r = r_new;
    end
end