%% Math 135A Project 1 - Joseph Benites jbeni047@ucr.edu

% -----------------------------------------------------------------
% Problem 1: Solve the polynomial Taylor Series using Horners Algorithm
% Write down the Taylor series based on the obtained coefficients.
% Given: p(x) = x^6 + 8x^5 - 3x^3 - 4x^2 + 10, about -3.
% -----------------------------------------------------------------

fprintf('---Problem 1 Results ---');
fprintf('\n\n');

%Define the polynomial
Poly = [ 1, 8, 0, -3, -4, 0, 10];

%Divisor define here for (x - (-3) = (x + 3)
Div = [1, 3];

%Define to store the coefficients
Coe = [];

%Loop Iteration here
for i = 1:length(Poly)
%Quotient Remainder and polynomial division
    [Quo, Rem] = deconv(Poly, Div);

    %Store remainder in Coe
    Coe = [Coe, Rem(end)];
    %New polynomial next iteration
    Poly = Quo;
end


disp('Taylor Coefficients (c_0 to c_6):');
disp(Coe);

fprintf('\nFinal Taylor Series about a = -3:\n');
fprintf('p(x) = -1160 + 1725(x+3) - 922(x+3)^2 + 177(x+3)^3 + 15(x+3)^4 - 10(x+3)^5 + 1(x+3)^6\n');


% -----------------------------------------------------------------
%Problem 2: Decimal to Octal conversion using standard division and
%multiplication
%Given: 42871.02044249548245602455 convert to octal into a single octal
%number/string
% -----------------------------------------------------------------

fprintf('---Problem 2 Results ---');
fprintf('\n\n');

% Part A: Integer section
int_part = 42871;
octal_int_digits = []; %To store our results

while int_part > 0
    %Store initial remainder
    digit = rem(int_part, 8);
        
    %Store into octal digit and in correct order.
    octal_int_digits = [digit, octal_int_digits];

    %New int
    int_part = floor(int_part / 8);
end

%Part B: Decimal Section

frac_part = 0.02044249548245602455;
octal_fract_list = [];
max_iteration = 20; %just so it does not go on for too long.
count = 0;

while frac_part > 0 && count < max_iteration
    % Multiply the fractional part by 8
    frac_part = frac_part * 8;
    
    % Store the integer part of the result
    digit = floor(frac_part);
    octal_fract_list = [octal_fract_list, digit];
    
    % Update the fractional part
    frac_part = frac_part - digit;
    
    % Increment the count
    count = count + 1;
end

disp('Integer Part in Octal: ');
disp(octal_int_digits);

disp('Decimal Part in Octal: ');
disp(octal_fract_list);

% Combine integer and fractional parts into a single octal string
octal_result = strcat(num2str(octal_int_digits), '.', num2str(octal_fract_list));
fprintf('\nOctal representation: \n%s\n', octal_result);



