clc;

A = input('Enter your coefficient Matrix: ');  % Coefficient matrix
B = input('Enter source vector: ');   % Right-hand side vector

% Ensure B is a column vector
if size(B, 1) ~= size(A, 1)
    error('The length of the source vector must match the number of rows in the coefficient matrix.');
end

N = length(B);
X = zeros(N, 1);

aug = [A B];  % Augmented matrix

operation_count_total = 0;  % Initialize counter for total number of operations
operation_count_upper = 0;  % Initialize counter for operations in reduction to upper triangular form

% Forward elimination with partial pivoting
for j = 1:N-1
    % Perform partial pivoting
    [~, max_idx] = max(abs(aug(j:N, j)));  % Find the row with the largest absolute value in the current column
    max_idx = max_idx + j - 1;  % Adjust index relative to the full matrix
    
    % Display the max_idx for each operation
    fprintf('For column %d, the maximum index is %d\n', j, max_idx);
    
    if max_idx ~= j
        % Swap the current row with the row containing the largest element in the column
        aug([j, max_idx], :) = aug([max_idx, j], :);
        
        % Display the row swap operation
        fprintf('Swapped row %d with row %d for partial pivoting\n', j, max_idx);
        disp('Updated Augmented Matrix after swapping:');
        disp(aug);
    end
    
    % Continue with elimination
    for i = j+1:N
        m = aug(i, j) / aug(j, j);  % Multiplier
        aug(i, :) = aug(i, :) - m * aug(j, :);  % Row operation
        
        % Count operations for upper triangular reduction
        operation_count_upper = operation_count_upper + 1;  % for the division (m = aug(i,j) / aug(j,j))
        operation_count_upper = operation_count_upper + 1;  % for the multiplication (m * aug(j,:))
        operation_count_upper = operation_count_upper + N;  % for the N subtractions (aug(i,:) - m * aug(j,:))
        
        % Display the operation details
        fprintf('Operation %d: Row %d = Row %d - (%.2f) * Row %d\n', ...
            operation_count_upper, i, i, m, j);
        
        % Display the updated augmented matrix after each row operation
        disp('Updated Augmented Matrix:');
        disp(aug);
    end
end

% Backward substitution
X(N) = aug(N, N + 1) / aug(N, N);  % Last variable
operation_count_total = operation_count_total + 1;  % Count the division for last variable

for k = N-1:-1:1
    % Count operations for each variable being calculated
    operation_count_total = operation_count_total + (N - k);  % Multiplications for (aug(k, k+1:N) * X(k+1:N))
    X(k) = (aug(k, N + 1) - aug(k, k + 1:N) * X(k + 1:N)) / aug(k, k);  % Corrected order of operations
    operation_count_total = operation_count_total + 1;  % Count the division for X(k)
end

% Display the solution vector X
disp('Solution X:');
disp(X);

% Display the total number of floating-point operations for upper triangular reduction
disp('Total number of floating-point operations for upper triangular reduction:');
disp(operation_count_upper);

% Display the total number of floating-point operations performed
disp('Total number of floating-point operations performed:');
disp(operation_count_total + operation_count_upper);  % Sum of upper triangular reduction and backward substitution
