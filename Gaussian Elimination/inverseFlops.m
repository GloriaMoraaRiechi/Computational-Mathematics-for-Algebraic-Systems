clc;

A = input('Enter the matrix to find its inverse: ');

A = sym(A);

% Check if the matrix is square
[n, m] = size(A);
if n ~= m
    error('Matrix must be square to find its inverse.');
end

% Augment A with the identity matrix
aug = [A eye(n)];

% Initialize FLOP counter
flops = 0;

% Perform Gaussian elimination (Gauss-Jordan elimination)
for j = 1:n
    % Find the pivot element and check if it's non-zero
    if aug(j,j) == 0
        % Perform partial pivoting if the pivot element is zero
        [~, max_idx] = max(abs(aug(j:n, j)));
        max_idx = max_idx + j - 1;  % Adjust index relative to full matrix
        if max_idx ~= j
            aug([j, max_idx], :) = aug([max_idx, j], :);  % Swap rows
            fprintf('Swapped row %d with row %d for partial pivoting\n', j, max_idx);
        end
    end
    
    % Normalize the pivot row (make pivot element equal to 1)
    pivot = aug(j,j);
    aug(j, :) = aug(j, :) / pivot;
    flops = flops + 2 * n;  % n divisions + n multiplications for each row
    
    % Eliminate other rows to make the column elements above and below the pivot zero
    for i = 1:n
        if i ~= j
            multiplier = aug(i,j);
            aug(i, :) = aug(i, :) - multiplier * aug(j, :);
            flops = flops + 2 * n;  % n multiplications + n subtractions
        end
    end
end

% Extract the inverse from the augmented matrix
A_inv = aug(:, n+1:end);

% Display the inverse matrix
disp('The inverse matrix is:');
disp(A_inv);

% Display the calculated FLOPs
fprintf('The total number of FLOPs for matrix inversion is: %.0f FLOPs\n', flops);
