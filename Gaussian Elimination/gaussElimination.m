clc

A = input('Enter the coefficient matrix: ');
B = input('Enter the source vector: ');

N = length(B);
X = zeros(N, 1);
Aug = [A B];  % Augmented matrix

% Forward elimination with partial pivoting
for j = 1:N-1
    % Find the row with the largest element in column j
    [~, maximumIndex] = max(abs(Aug(j:N, j)));
    maximumIndex = maximumIndex + j - 1;  % Adjust to full matrix index

    % Display the maximum index
    fprintf('For column %d, maximum index is %d\n', j, maximumIndex);

    % Row swap if necessary
    if maximumIndex ~= j
        Aug([j, maximumIndex], :) = Aug([maximumIndex, j], :);  % Swap rows
        disp(Aug);  % Display the augmented matrix after each swap
    end
    
    % Forward elimination (Gaussian elimination)
    for i = j+1:N
        m = Aug(i, j) / Aug(j, j);  % Multiplier
        Aug(i, :) = Aug(i, :) - m * Aug(j, :);  % Row operation
    end
end

% Backward substitution
X(N) = Aug(N, N+1) / Aug(N, N);  % Solve for the last variable
for k = N-1:-1:1
    X(k) = (Aug(k, N+1) - Aug(k, k+1:N) * X(k+1:N)) / Aug(k, k);  % Solve for X(k)
end

disp(X);
