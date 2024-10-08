clc;

A = input('Enter the matrix to find its inverse: ');

% Check if the matrix is square
[n, m] = size(A);
if n ~= m
    error('Matrix must be square to find its inverse.');
end

% Augment A with the identity matrix
aug = [A eye(n)];

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
    aug(j, :) = aug(j, :) / aug(j,j);
    
    % Eliminate other rows to make the column elements above and below the pivot zero
    for i = 1:n
        if i ~= j
            m = aug(i,j);
            aug(i, :) = aug(i, :) - m * aug(j, :);
        end
    end
end

% Extract the inverse from the augmented matrix
A_inv = aug(:, n+1:end);

% Display the inverse matrix
disp('The inverse matrix is:');
disp(A_inv);
