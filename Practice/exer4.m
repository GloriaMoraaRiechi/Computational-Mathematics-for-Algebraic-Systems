% Define the matrix A (for example purposes, this will be randomly generated)
% Replace this with your actual matrix A
A = randi(100, 8, 9);  % Example 8x9 matrix with random integers from 1 to 100

% Display the original matrix A
disp('Original matrix A:');
disp(A);

% Initialize variables
[m, n] = size(A);          % Get the size of matrix A
resultMatrix = zeros(m, n); % Preallocate matrix to maximum possible size
resultCount = 0;           % Initialize counter for the number of rows added

% Iterate over each row of matrix A
for i = 1:m
    row = A(i, :);         % Extract the i-th row
    hasMultipleOf5 = any(mod(row, 5) == 0);  % Check for multiples of 5
    hasMultipleOf7 = any(mod(row, 7) == 0);  % Check for multiples of 7
    
    % If the row contains both multiples of 5 and 7, add it to resultMatrix
    if hasMultipleOf5 && hasMultipleOf7
        resultCount = resultCount + 1; % Increment the counter
        resultMatrix(resultCount, :) = row; % Store the row
    end
end

% Remove any preallocated but unused rows
if resultCount > 0
    resultMatrix = resultMatrix(1:resultCount, :);
else
    resultMatrix = [];  % If no rows meet the criteria, return an empty matrix
end

% Display the resulting matrix
disp('The matrix containing rows with at least one multiple of 5 and one multiple of 7 is:');
disp(resultMatrix);
