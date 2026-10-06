%% --- MATH135A Project 2: Problem 1 ---
% Solves a 4x4 system using Gaussian elimination with scaled partial pivoting


% --- 1. Setup ---
A = [ 0.4096, 0.1234, 0.3678, 0.2943;
      0.2246, 0.3872, 0.4015, 0.1129;
      0.3645, 0.1920, 0.3781, 0.0643;
      0.1784, 0.4002, 0.2786, 0.3927];

B = [ 0.4042; % Using the value 0.4042 as you provided
      0.1250;
      0.4255;
      0.2557];

% Get the size of the matrix
n = 4;

% Create the augmented matrix
Aug = [A, B];

% Create the scale vector (as a column)
% We can also calculate this in MATLAB:
s = max(abs(A), [], 2); % Finds max absolute value in each row

disp('Initial Augmented Matrix:');
disp(Aug);
disp('Initial Scale Vector:');
disp(s);

% --- 2. Forward Elimination ---
for k = 1:n-1
    
    % --- Pivoting ---
    % Find scaled ratios for rows k through n in column k
    ratios = abs(Aug(k:n, k)) ./ s(k:n);
    
    % Find the index of the largest ratio
    [max_val, max_idx] = max(ratios);
    
    % Calculate the actual row number to pivot with
    pivot_row = k + max_idx - 1;
    
    % Swap rows if pivot_row is different from the current row k
    if pivot_row ~= k
        % Swap rows in Aug matrix
        temp_row = Aug(k, :);
        Aug(k, :) = Aug(pivot_row, :);
        Aug(pivot_row, :) = temp_row;
        
        % IMPORTANT: Swap elements in the scale vector too
        temp_s = s(k);
        s(k) = s(pivot_row);
        s(pivot_row) = temp_s;
    end

    % --- Elimination ---
    % Loop through all rows below the pivot row k
    for i = k+1:n
        
        % Find the multiplier
        m = Aug(i, k) / Aug(k, k);
        
        % Subtract the scaled pivot row from the current row i
        Aug(i, :) = Aug(i, :) - m * Aug(k, :);
    end
end

disp('Upper Triangular (Eliminated) Matrix:');
disp(Aug);

% --- 3. Back Substitution ---
% Create an empty solution vector (filled with zeros)
x = zeros(n, 1); 

% Solve for the last element, x(n)
x(n) = Aug(n, n+1) / Aug(n, n);

% Loop backwards from n-1 down to 1
for k = n-1:-1:1
    
    % Calculate the sum of the parts we already know
    % (Aug(k, k+1:n) is a row vector, x(k+1:n) is a column vector)
    sum_of_known_parts = Aug(k, k+1:n) * x(k+1:n);
    
    % Solve for x(k)
    x(k) = (Aug(k, n+1) - sum_of_known_parts) / Aug(k, k);
    
end

% --- 4. Final Answer ---
disp('Final Answer Vector x:');
disp(x);