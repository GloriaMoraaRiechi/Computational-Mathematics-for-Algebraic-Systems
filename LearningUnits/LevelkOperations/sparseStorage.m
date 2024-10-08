clc
% Define the matrix A
A = [1 0 3 0; 0 1 0 4];

% Convert A to sparse form
A_sparse = sparse(A);

% Display the sparse matrix
disp(A_sparse);

[row, col, val] = find(A_sparse);
disp('Row Indices:');
disp(row);
disp('Column Indices:');
disp(col);
disp('Values:');
disp(val);

ASparseT = sparse(A');
disp(ASparseT)

%% 
clc
A = cell(2,1) %create a cell array with 2 rows
% Define sparse matrix A in MATLAB using 1-based indexing
A{1} = [1, 1; 4, 2]  % Row 1: A(1,1) = 1, A(1,4) = 2
A{2} = [2, 3; 3, 4]  % Row 2: A(2,2) = 3, A(2,3) = 4

% Display the sparse matrix A
disp('Sparse Matrix A:');
for i = 1:numel(A)
    disp(A{i});
end


%% 
% propose a sparse matrix multiplication algorithm that is based on the row combination interpretation. assume both matrices are stored by rows

clc
A = cell(2,1);
B = cell(2,1);
function C = sparse_matrix_multiply(A, B)
    % A and B are cell arrays representing sparse matrices stored by rows.
    
    m_A = numel(A);  % Number of rows in A
    n_B = max(cellfun(@(x) max(x(:,1)), B));  % Number of columns in B
    C = zeros(m_A, n_B);  % Initialize the result matrix
    
    % Iterate through each row of A
    for i = 1:m_A
        % For each non-zero entry in row i of A
        for j = 1:size(A{i}, 1)
            j_A = A{i}(j, 1);  % Column index in A
            val_A = A{i}(j, 2);  % Value in A
            
            % Check if j_A is a valid index for B
            if j_A <= numel(B)
                % Iterate through the corresponding row in B
                for k = 1:size(B{j_A}, 1)
                    j_B = B{j_A}(k, 1);  % Column index in B
                    val_B = B{j_A}(k, 2);  % Value in B
                    
                    % Update the result matrix C
                    C(i, j_B) = C(i, j_B) + val_A * val_B;  % Accumulate the product
                end
            else
                fprintf('Warning: Index %d for B is out of bounds\n', j_A);
            end
        end
    end
    
    % Display the result matrix
    disp('Result of Ab:');
    disp(C);
end

% Example usage
A{1} = [1, 1; 4, 2];  % Row 1: A(1,1) = 1, A(1,4) = 2
A{2} = [2, 3; 3, 4];  % Row 2: A(2,2) = 3, A(2,3) = 4

B{1} = [1, 4; 3, 5];  % Row 1: B(1,1) = 4, B(1,3) = 5
B{2} = [2, 1; 4, 6];  % Row 2: B(2,2) = 1, B(2,4) = 6

% Perform the sparse matrix multiplication
C = sparse_matrix_multiply(A, B);
