%% --- MATH135A Project 2: Problem 2 ---
% Solves a tridiagonal system for n=100 using Procedure Tri



disp('--- Problem 2 (n=100) ---');
n2 = 100;

% Setup the RHS vector 'b'
b = 40 * ones(n2, 1);
b(1) = -20;
b(end) = -20;

% Setup the three diagonals
% We make them size n, though a(1) and c(n) won't be used
a = -1 * ones(n2, 1); % Sub-diagonal
d =  4 * ones(n2, 1); % Main diagonal
c = -1 * ones(n2, 1); % Super-diagonal

% Pre-allocate solution vector
x2 = zeros(n2, 1);

% --- 2. Forward Elimination (Thomas Algorithm) ---
for i = 2:n2
    % Find the multiplier
    m = a(i) / d(i-1);
    
    % Update the main diagonal (d) and the RHS (b)
    d(i) = d(i) - m * c(i-1);
    b(i) = b(i) - m * b(i-1);
end

% --- 3. Back Substitution (Thomas Algorithm) ---
% Solve for the last x
x2(n2) = b(n2) / d(n2);

% Loop backwards from n-1 down to 1
for k = n2-1:-1:1
    % Solve for x(k)
    x2(k) = (b(k) - c(k) * x2(k+1)) / d(k);
end

% --- 4. Final Answer (Problem 2) ---
disp('Final Answer Vector x for Problem 2 (n=100):');
disp(x2);